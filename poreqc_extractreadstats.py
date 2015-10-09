#!/usr/bin/env python

# ============================================================================ #
# poreqc_extractreadstatsfrompostfiles.py                                      #
# ============================================================================ #
"""
Extract MinION event read statistics

This program just extracts the event data from the basecalled fast5 files.
Really, there should now only be one program to extract event and called
stats, but that will take too long to re-program so I'm just replacing this
part of the pipeline.
"""

# ============================================================================ #
# Import Modules                                                               #
# ============================================================================ #

import argparse, datetime, h5py, math, numpy as np, os, re, shlex, stat, sys, time

# ============================================================================ #
# Global variables                                                             #
# ============================================================================ #

_progdesc = 'Extract MinION event read statistics'

_ErrorInvalidData = 18
_ErrorOutfileExists = 31
_ErrorFileRename = 33

_rawread_dtype = [
    # From filename _chN_fileN.fast5
    ('fast5_filename', 'S100'),
    ('run_number', 'S100'),
    ('file_number', np.int),
    # Key/read_id
    ('channel_number', np.int),
    ('read_number', np.int),
    # Key/tracking_id
    ('asic_id', 'S100'),
    ('asic_temp', np.float),
    ('device_id', 'S100'),
    ('exp_script_purpose', 'S100'),
    ('exp_start_time', np.int),
    ('exp_start_localtime_iso', 'S100'),
    ('flow_cell_id', 'S100'),
    ('heatsink_temp', np.float),
    ('run_id', 'S100'),
    ('version_name', 'S100'),
    # Reads/Read_READ_NUMBER
    ('event_count', np.int),
    ('median_before', np.float),
    ('sample_rate', np.float),
    ('start', np.int),
    ('start_mux', np.int),
    # Sequences/Meta
    ('numerical_encoding', 'S20'),
    ('precision', 'S20'),
    ('tool', 'S100'),
    ('version', 'S20'),
    ('read_start_seconds', np.float),
    ('read_end_seconds', np.float),
    ('read_start_localtime_iso', 'S100'),
    ('read_end_localtime_iso', 'S100'),
    ('read_start_hours_since_expt_start', np.float),
    ('read_start_hours_since_expt_start_qtrhrbucketname', np.float),
    ('read_end_hours_since_expt_start', np.float),
    ('read_end_hours_since_expt_start_qtrhrbucketname', np.float),
    ('read_end_hours_since_read_start', np.float),
    ('read_end_hours_since_read_start_qtrhrbucketname', np.float),
    ('read_duration_seconds', np.float),
    ('events_per_second', np.float),
    ('eps_1sthalf', np.float),
    ('eps_2ndhalf', np.float)
]
_read_dtype = None	# The dtype with all the eps_b1-10 values as well
_read_dtype = _rawread_dtype[:]
for i in range(0, 10):
    _read_dtype.append(('eps_b{0}'.format(i+1), np.float))

_readstatsfile_was_updated = False

# ============================================================================ #
# Program usage                                                                #
# ============================================================================ #

def Initialise():
    'Process command-line arguments and define output filenames.'

    global _progdir, _progname, _args
    _progdir = os.path.dirname(os.path.realpath(sys.argv[0]))
    _progname = os.path.basename(os.path.realpath(sys.argv[0]))
    progexamples = []
    parser = argparse.ArgumentParser(description=_progdesc, \
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--inreadsdir', metavar='str', dest='inreadsdir', \
        type=str, default=None, help='Input directory containing the raw .fast5 files of event data', required=True)
    parser.add_argument('--inreadsarecalled', action='store_true', dest='inreadsarecalled', \
        help='Set to true if inreadsdir contains basecalled fast5 files', required=False)
    parser.add_argument('--outstatdir', metavar='str', dest='outstatdir', \
        type=str, default=None, help='Output directory for .txt files of read statistics', required=True)
    parser.add_argument('--outprefix', metavar='str', dest='outprefix', \
        type=str, default=None, help='All output files will start with this prefix', required=True)
    parser.add_argument('--overwrite', action='store_true', dest='overwrite', \
        help='Force recomputation of all output files, overwriting existing as necessary', required=False)
    parser.add_argument('--debug', action='store_true', dest='debug', \
        help='Print verbose diagnostics (value must be 0 or 1)', required=False)
    _args = parser.parse_args()

    global _single_path,  _simplified_path, _filldelay_path, _activechannels_path
    _single_path = os.path.join(_args.outstatdir, _args.outprefix + '_readstats.txt')
    _simplified_path = os.path.join(_args.outstatdir, _args.outprefix + '_readstats_simplified.txt')
    _filldelay_path = os.path.join(_args.outstatdir, _args.outprefix + '_filldelay.txt')
    _activechannels_path = os.path.join(_args.outstatdir, _args.outprefix + '_activechannels.txt')

def Get_RawFast5Version(hdf):
    'Return a string describing the raw fast5 version of the file.'
    version = 'unknown'
    if len(hdf.keys()) == 3 and 'Key' in hdf.keys() and 'Reads' in hdf.keys() and 'Sequences' in hdf.keys():
        version = 'raw_v1'
    elif len(hdf.keys()) == 3 and 'Analyses' in hdf.keys() and 'Sequences' in hdf.keys() and 'UniqueGlobalKey' in hdf.keys():
        version= 'raw_v2'
    return version

def Get_ReadStatsFromFile(fast5_path):
    'Parse the fast5 events file and return a np.array of the statistics.'
    fast5_filename = os.path.basename(fast5_path)
    if os.stat(fast5_path).st_size == 0:
        sys.stdout.write('Warn: Ignoring empty file ({0})\n'.format(fast5_path))
        return None
    run_number = '_'.join(fast5_filename.split('_')[-5:-3])
    file_number = int(fast5_filename.split('_')[-2].replace('file', ''))
    try:
        hdf = h5py.File(fast5_path, 'r')
    except:
        sys.stdout.write('Warn: Ignoring file that could not be opened with h5py.File ({0})\n'.format(fast5_path))
        return None

    # The data format changes - both the 'directory' structure and where attributes live.
    # So need to add some logic here to work out where to find things.
    raw_fast5_version = Get_RawFast5Version(hdf)

    if raw_fast5_version == 'raw_v1':
        channel_number =  int(hdf['Key/read_id'].attrs['channel_number'])
        read_number =  int(hdf['Key/read_id'].attrs['read_number'])
        key2 = 'Key/tracking_id'
        key3 = 'Reads/Read_{0}'.format(read_number)
        key4 = 'Sequences/Meta'
        keyevents = 'Reads/Read_{0}/Events'.format(read_number)
        startidx = 0 if not _args.inreadsarecalled else 2
        # CI 2014-08-18: This fails sometimes because the hdf[keyevents][()] returns error
        # "OverflowError: cannot fit 'long' into an index-sized integer". I don't have a
        # fix for this so the solution is simply to discard this read.
        try:
            event_count = min(int(hdf[key3].attrs['event_count']), len(hdf[keyevents][()]))
        except:
            sys.stdout.write('Warn: Failed to retrieve event_count, probably due to bug in workflow software ({0})\n'.format(fast5_path))
            return None
        sample_rate = float(hdf[key3].attrs['sample_rate'])
        start = hdf[key3].attrs['start']
        start_mux = hdf[key3].attrs['start_mux']

    elif raw_fast5_version == 'raw_v2':
        channel_number =  int(hdf['UniqueGlobalKey/channel_id'].attrs['channel_number'])
        read_number = int(hdf['Analyses/EventDetection_000/Reads'].keys()[0].split('_')[1])
        key2 = 'UniqueGlobalKey/tracking_id'
        key3 = 'Analyses/EventDetection_000/Reads/Read_{0}'.format(read_number)
        key4 = 'Sequences/Meta'
        keyevents = 'Analyses/EventDetection_000/Reads/Read_{0}/Events'.format(read_number)
        startidx = 0 if not _args.inreadsarecalled else 2
        event_count = hdf['Analyses/EventDetection_000/Reads/Read_{0}/Events'.format(read_number)].len()
        sample_rate = float(hdf['UniqueGlobalKey/channel_id'].attrs['sampling_rate'])
        start = hdf[key3].attrs['start_time']
        start_mux = hdf[key3].attrs['start_mux']

    elif raw_fast5_version == 'unknown':
        sys.stdout.write('Erro: Unrecognised raw fast5 file format - please investigate\n')

    try:
        exp_start_time = float(hdf[key2].attrs['exp_start_time'])
        exp_start_time_localtime_iso = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(exp_start_time))
    except:
        sys.stdout.write('Warn: Failed to extract exp_start_time ({0})\n'.format(fast5_path))
        return None

    try:
        # CI 2015-07-09: In pre-basecalled, start is in column 1 (start, length, mean, variance),
        # in post-basecalled in column 3 (mean, stdv, start, length)
        read_start_seconds = exp_start_time + float(hdf[keyevents][()][0][startidx]) / sample_rate
        read_end_seconds = exp_start_time + float(hdf[keyevents][()][-1][startidx]) / sample_rate

        read_start_localtime_iso = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(read_start_seconds))
        read_end_localtime_iso = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(read_end_seconds))

        read_start_hours_since_expt_start = hdf[keyevents][()][0][startidx] / sample_rate / 3600.0
        read_start_hours_since_expt_start_qtrhrbucketname = '{0:.2f}'.format(math.ceil(read_start_hours_since_expt_start / 0.25) * 0.25)

        read_end_hours_since_expt_start = hdf[keyevents][()][-1][startidx] / sample_rate / 3600.0
        read_end_hours_since_expt_start_qtrhrbucketname = '{0:.2f}'.format(math.ceil(read_end_hours_since_expt_start / 0.25) * 0.25)

        read_end_hours_since_read_start = (hdf[keyevents][()][-1][startidx] - hdf[keyevents][()][0][startidx]) / sample_rate / 3600.0
        read_end_hours_since_read_start_qtrhrbucketname = '{0:.2f}'.format(math.ceil(read_end_hours_since_read_start / 0.25) * 0.25)

        read_duration_seconds = read_end_seconds - read_start_seconds
    except:
        sys.stdout.write('Warn: Failed to read_start-/read_end_seconds, probably due to bug in workflow software producing an event array that cannot be accessed ({0})\n'.format(fast5_path))
        return None

    # Example of file that has event_count != len(events)
    # BOWDEN03_MdC_S.aureus_run_2.1_BN_011_34887GT_5913_1_ch101_file0_strand.fast5
    #event_count = min(int(hdf[key3].attrs['event_count']), len(hdf[keyevents][()]))
    if event_count > 100:
        events_per_second = (event_count / float(read_duration_seconds)) if (read_duration_seconds) else 0.0
        events_in_1sthalf = int(event_count / 2.0)
        num_seconds_in_1sthalf = (float(hdf[keyevents][()][events_in_1sthalf][startidx]) - float(hdf[keyevents][()][0][startidx])) / sample_rate
        num_seconds_in_2ndhalf = (float(hdf[keyevents][()][-1][startidx]) - float(hdf[keyevents][()][events_in_1sthalf+1][startidx])) / sample_rate
        eps_1sthalf = (events_in_1sthalf / num_seconds_in_1sthalf) if num_seconds_in_1sthalf else 0.0
        eps_2ndhalf = ((event_count - events_in_1sthalf) / num_seconds_in_2ndhalf) if num_seconds_in_2ndhalf else 1.0
    
        eps_boundaries = np.arange(0, event_count, int(round(event_count/10.0)))
        if len(eps_boundaries) == 10:
            eps_boundaries = np.append(eps_boundaries, event_count - 1)
        elif len(eps_boundaries) == 11:
            eps_boundaries[-1] = event_count - 1
        else:
            sys.stderr.write('Erro: computed {0} eps_boundaries instead of exactly 11 entries in array - investigate.\n'.format(len(eps_boundaries)))
            sys.exit(_ErrorInvalidData)
        eps_bucket = [0.0] * (len(eps_boundaries)-1)
        for i in range(1, len(eps_boundaries.tolist())):
            num_events = eps_boundaries[i] - eps_boundaries[i-1]
            try:
                num_seconds = float(float(hdf[keyevents][()][eps_boundaries[i]][startidx]) - float(hdf[keyevents][()][eps_boundaries[i-1]][startidx])) / sample_rate
                eps_bucket[i-1] = num_events / num_seconds
            except:
                eps_bucket[i-1] = 0.0
            pass
    else:
        events_per_second = 0.0
        events_in_1sthalf = 0
        num_seconds_in_1sthalf = 0
        num_seconds_in_2ndhalf = 0
        eps_1sthalf = 0.0
        eps_2ndhalf = 1.0
        eps_bucket = [0.0] * 10

    data_fromfast5 = [
        fast5_filename,
        run_number,
        file_number,
        channel_number,
        read_number,
        hdf[key2].attrs['asic_id'],
        float(hdf[key2].attrs['asic_temp']),
        hdf[key2].attrs['device_id'],
        hdf[key2].attrs['exp_script_purpose'],
        int(hdf[key2].attrs['exp_start_time']),
        exp_start_time_localtime_iso,
        hdf[key2].attrs['flow_cell_id'],
        float(hdf[key2].attrs['heatsink_temp']),
        hdf[key2].attrs['run_id'],
        hdf[key2].attrs['version_name'],
        event_count,
        hdf[key3].attrs['median_before'],
        sample_rate,
        start,
        start_mux,
        ''.join(hdf[key4].attrs['numerical_encoding']),
        ''.join(hdf[key4].attrs['precision']),
        ''.join(hdf[key4].attrs['tool']),
        ''.join(hdf[key4].attrs['version']),
        read_start_seconds,
        read_end_seconds,
        read_start_localtime_iso,
        read_end_localtime_iso,
        read_start_hours_since_expt_start,
        float(read_start_hours_since_expt_start_qtrhrbucketname),
        read_end_hours_since_expt_start,
        float(read_end_hours_since_expt_start_qtrhrbucketname),
        read_end_hours_since_read_start,
        float(read_end_hours_since_read_start_qtrhrbucketname),
        read_duration_seconds,
        events_per_second,
        eps_1sthalf,
        eps_2ndhalf
    ]
    hdf.close()
    rowL = [str(x) for x in data_fromfast5 + eps_bucket]

    return rowL

def Readstats_FileToArray():
    'Read entire readstats.txt file into globals _E and _Efile.'
    global _E, _Efile
    _E = np.loadtxt(_single_path, dtype=_read_dtype, delimiter='\t', skiprows=1)
    _Efile = {}
    for i in range(0, len(_E)):
        _Efile[_E[i][0]] = 1

def ReadStats_Fast5ToFile(indir, filelist):
    'Write all existing data from _E to readstats.txt, then append rows from new files in filelist.'

    # If readstats.txt file already exists, rename it so it is can be accessed later if necessary.'
    if os.path.exists(_single_path):
        new_single_path = _single_path.replace('.txt', '_' + time.strftime('%Y%m%d-%H%M%S', time.localtime()) + '.txt')
        try:
            os.rename(_single_path, new_single_path)
        except:
            sys.stderr.write('Error: Failed to make backup existing readstats.txt file ({0} to {1})\n'.format(_single_path, new_single_path))
            sys.exit(_ErrorFileRename)

    # Open the output readstats.txt file.
    lenfilelist = len(filelist)
    single_fp = open(_single_path, 'w')
    single_fp.write('{0}\n'.format('\t'.join([x[0] for x in _read_dtype])))
    # dump all the existing data from _E into it.
    if _E is not None:
        for i in range(0, len(_E)):
            row = [str(x) for x in _E[i]]
            single_fp.write('{0}\n'.format('\t'.join(row)))
        single_fp.flush()
    for j, file in enumerate(filelist):
        if _args.debug and (j+1) >= 500:
            break   # XXXX DEBUGGING
        thing_path = os.path.join(indir, file)
        if _Efile.has_key(file):
            continue
        row = Get_ReadStatsFromFile(thing_path)
        if row is not None:
            single_fp.write('{0}\n'.format('\t'.join(row)))
            single_fp.flush()
            if (((j+1) % 100) == 0) or ((j+1) == lenfilelist):
                sys.stdout.write('{0} Info: Processing {1}-th of {2} additional fast5 file\n'.format(
                    time.strftime('%Y%m%d-%H%M%S', time.localtime()), j+1, lenfilelist))
                single_fp.flush()
        else:
            sys.stdout.write('{0} Info: Processing failed for {1}-th additional fast5 file ({2})\n'.format(
                time.strftime('%Y%m%d-%H%M%S', time.localtime()), j+1, thing_path))
    # And close the readstats.txt file.
    if single_fp:
        single_fp.close()

def Extract_ReadStats():
    '''
    Iterate through all reads/prefix.fast5 files and output outdir/prefix_readstats.txt files (if necessary).
    Pre-conditions:
        - There might be an existing _single_path but it might be incomplete.
    Post-conditions:
        - readstats.txt file and _E array contain all the available information from the raw fast5 files.
    '''

    # Set up IO directories.
    inreadsdir = os.path.expandvars(_args.inreadsdir)
    outstatdir = os.path.expandvars(_args.outstatdir)

    # Initialise global variables.
    global _E, _Efile
    _E = None
    _Efile = {}

    # Get current list of raw fast5 files to be processed.
    files = [x for x in os.listdir(inreadsdir) if os.path.isfile(os.path.join(inreadsdir, x)) and x.endswith('.fast5')]
    files.sort()

    # If readstats.txt file exists ...
    #     check if readstats.txt file is complete and go to next section if it is.
    # ... if readstats.txt does not exist, update with info from all raw fast5 files.
    if os.path.exists(_single_path):
        Readstats_FileToArray()
        numrecords_before = len(_E)
        if len(_E) == len(files) and not _args.overwrite:
            sys.stdout.write('Info: readstats.txt is already up-to-date, no need to compute ({0})\n'.format(_single_path))
            numrecords_after = numrecords_before
        else:
            if not _args.overwrite:
                sys.stdout.write('Info: readstats.txt exists but may not be complete, starting update now ({0})\n'.format(_single_path))
            else:
               _E = None
               _Efile = {}
            ReadStats_Fast5ToFile(inreadsdir, files)
            Readstats_FileToArray()
            numrecords_after = len(_E)
    else:
        numrecords_before = 0
        sys.stdout.write('Info: readstats.txt does not exist, starting iteration over {0} fast5 files now ({1})\n'.format(len(files), _single_path))
        ReadStats_Fast5ToFile(inreadsdir, files)
        Readstats_FileToArray()
        numrecords_after = len(_E)
    global _readstatsfile_was_updated
    _readstatsfile_was_updated = (numrecords_before < numrecords_after) or _args.overwrite
    sys.stdout.write('{0} Info: readstats.txt updated with {1} new records\n'.format(time.strftime('%Y%m%d-%H%M%S', time.localtime()), numrecords_after - numrecords_before))

    # Create new readstats-simplified.txt file from all _E records, if required.
    if _readstatsfile_was_updated:
        if os.path.exists(_simplified_path):
            sys.stdout.write('{0} Info: Overwriting simplified file ({1})\n'.format(time.strftime('%Y%m%d-%H%M%S', time.localtime()), _simplified_path))
        else:
            sys.stdout.write('{0} Info: Creating new simplified file ({1})\n'.format(time.strftime('%Y%m%d-%H%M%S', time.localtime()), _simplified_path))
        with open(_simplified_path, 'w') as simplified_fp:
            simplified_fp.write('{0}\n'.format('\t'.join([_read_dtype[0][0], str(_read_dtype[17][0])])))
            for i in range(0, len(_E)):
                simplified_fp.write('{0}\t{1}\n'.format(_E[i][0], _E[i][17]))
    else:
        sys.stdout.write('{0} Info: readstats_simplified.txt already up to date\n'.format(time.strftime('%Y%m%d-%H%M%S', time.localtime())))

def Extract_ActiveChannels():
    '''
    Read in the summary readstats.txt file and produce a file containing the time_bucket and num_active_channels.
    Re-compute the information for the file. If nothing needs to be done, don't update the file.
    '''

    global _activechannels_path
    if not _readstatsfile_was_updated and os.path.exists(_activechannels_path) and not _args.overwrite:
        sys.stdout.write('{0} Info: activechannels.txt already up to date\n'.format(time.strftime('%Y%m%d-%H%M%S', time.localtime())))
        return
    if os.path.exists(_activechannels_path):
        sys.stdout.write('{0} Info: Overwriting activechannels.txt ({1})\n'.format(time.strftime('%Y%m%d-%H%M%S', time.localtime()), _activechannels_path))
    else:
        sys.stdout.write('{0} Info: Creating new activechannels.txt ({1})\n'.format(time.strftime('%Y%m%d-%H%M%S', time.localtime()), _activechannels_path))

    # Compute the number of unique channels in each time bucket.
    interval = 0.25
    active_channels = {}
    data = _E
    for i in range(0, len(data)):
        if (0):	# Debugging
            print '{0}\t{1}\t{2}\t{3}'.format(
                data['channel_number'][i],
                data['read_start_hours_since_expt_start_qtrhrbucketname'][i],
                data['read_end_hours_since_expt_start_qtrhrbucketname'][i],
                str(np.arange(data['read_start_hours_since_expt_start_qtrhrbucketname'][i],
                    data['read_end_hours_since_expt_start_qtrhrbucketname'][i]+interval, interval)))
        for bucket in np.arange(data['read_start_hours_since_expt_start_qtrhrbucketname'][i], data['read_end_hours_since_expt_start_qtrhrbucketname'][i]+interval, interval):
            if not active_channels.has_key(bucket):
                active_channels[bucket] = []
            if data['channel_number'][i] not in active_channels[bucket]:
                active_channels[bucket].append(data['channel_number'][i])

    # Debugging
    if (0):
        for bucket in np.arange(0, max(active_channels.keys())+interval, interval):
            num_active = 0
            list_active = []
            if active_channels.has_key(bucket):
                num_active = len(active_channels[bucket])
                list_active = active_channels[bucket]
                list_active.sort()
            print '{0}\t{1}\t{2}'.format(bucket, num_active, list_active)

    # Save the number of active channels in each time-bucket to a file.
    H = ['time_bucket_since_expt_start_hrs', 'num_active_channels']
    with open(_activechannels_path, 'w') as active_fp:
        active_fp.write('{0}\n'.format('\t'.join(H)))
        for bucket in np.arange(0, max(active_channels.keys())+interval, interval):
            num_active = len(active_channels[bucket]) if active_channels.has_key(bucket) else 0
            R = [bucket, num_active]
            active_fp.write('{0}\n'.format('\t'.join([str(x) for x in R])))

def Extract_FillDelay():
    '''
    Create file with the fill_delay values (in seconds), or NA if the
    file_number within a channel is not contiguous.
    '''

    if not _readstatsfile_was_updated and os.path.exists(_filldelay_path) and not _args.overwrite:
        sys.stdout.write('{0} Info: filldelay.txt already up to date\n'.format(time.strftime('%Y%m%d-%H%M%S', time.localtime())))
        return
    if os.path.exists(_filldelay_path):
        sys.stdout.write('{0} Info: Overwriting filldelay.txt ({1})\n'.format(time.strftime('%Y%m%d-%H%M%S', time.localtime()), _filldelay_path))
    else:
        sys.stdout.write('{0} Info: Creating new filldelay.txt ({1})\n'.format(time.strftime('%Y%m%d-%H%M%S', time.localtime()), _filldelay_path))

    data = np.loadtxt(_single_path, dtype=_read_dtype, delimiter='\t', skiprows=1)
    L = []
    H = ['run_number', 'channel_number', 'file_number', 'read_number', 'sample_rate', 'exp_start_time', 'start', 'initfilldelay_seconds', 'refilldelay_seconds']
    for i in range(0, len(data)):
        row = (data[i][1], data[i][3], data[i][2], data[i][4], data[i][17], data[i][9], data[i][18], 0.0, 0.0)
        L.append(row)
    A = np.array(L, dtype=[('run_number', 'S20'), ('channel_number', np.int), ('file_number', np.int), ('read_number', np.int), ('sample_rate', np.int), ('exp_start_time', np.int), ('start', np.int), ('initfilldelay_seconds', np.float), ('refilldelay_seconds', np.float) ])
    A = np.sort(A, order=['run_number', 'channel_number', 'file_number', 'read_number'])
    prev_row = ['NA', 0, 0, 0, 0, 0, 0, 0.0, 0.0]
    prev_bn, prev_cn, prev_fn, prev_rn, prev_sr, prev_est, prev_s, prev_ifd, prev_rfd = prev_row
    
    # Insert the filldelay values, which are in the same units as start and length (i.e., 1/sample_rate of a second).
    for i in range(0, len(A)):
        row = A[i]
        bn, cn, fn, rn, sr, est, s, ifd, rfd = row
        if fn == 0:
            A[i][7] = s / float(sr)
        else:
            A[i][7] = -1.0
        if bn == prev_bn and cn == prev_cn and fn == (prev_fn+1) and fn > 0:
            A[i][8] = (s - prev_s) / sr
        else:
            A[i][8] = -1.0
        prev_row = row
        prev_bn, prev_cn, prev_fn, prev_rn, prev_sr, prev_est, prev_s, prev_ifd, prev_rfd = prev_row

    # Save the fill delay values to a file.
    result = list(A)
    with open(_filldelay_path, 'w') as filldelay_fp:
        filldelay_fp.write('{0}\n'.format('\t'.join(H)))
        for i in range(0, len(result)):
            row = list(result[i])
            if row[7] == -1:
                row[7] = 'NA'
            if row[8] == -1:
                row[8] = 'NA'
            filldelay_fp.write('{0}\n'.format('\t'.join([str(x) for x in row])))

    # Debugging
    if (0):
        print '{0}'.format('\t'.join(H))
        for i in range(0, len(A)):
            print '{0}'.format('\t'.join([str(x) for x in A]))

# ============================================================================ #
# Main                                                                         #
# ============================================================================ #

if __name__ == '__main__':
    Initialise()
    Extract_ReadStats()
    Extract_ActiveChannels()
    Extract_FillDelay()

# ============================================================================ #
