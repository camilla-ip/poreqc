#!/usr/bin/env python

# ============================================================================ #
# extract_callstats.py                                                         #
# ============================================================================ #
"""
Extract MinION base-call statistics
"""

# ============================================================================ #
# Import Modules                                                               #
# ============================================================================ #

import argparse, datetime, h5py, math, numpy as np, os, random, re, shlex, stat, sys, time

# ============================================================================ #
# Global variables                                                             #
# ============================================================================ #

_progdesc = 'Extract MinION base-call statistics'

_ErrorDirCreate = 13

# The simplified output from extract_readstats.py, containing information that
# should, really, be in the analysed basecall files but isn't.
_readstatssimp_dtype = [
    ('fast5_filename', 'S100'),
    ('sample_rate', np.float)
]

# This has to be generic enough to be applicable to the DNA (2D) and cDNA (1D) workflows.
_rawcall_dtype = [
    # From filename _chN_fileN.fast5
    ('fast5_filename', 'S100'),
    ('run_number', 'S100'),
    ('file_number', np.int),
    ('read_number', np.int),
    # UniqueGlobalKey/channel_id
    ('channel_number', np.int),
    ('sampling_rate', np.float),
    # UniqueGlobalKey/tracking_id
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
    # Sequences/Meta
    ('numerical_encoding', 'S20'),
    ('precision', 'S20'),
    ('tool', 'S100'),
    ('version', 'S20'),
    # EventDetection_000/Reads/Read_N
    ('duration', np.float),
    ('start_time', np.float),
    ('read_start_localtime_iso', 'S100'),
    ('read_end_localtime_iso', 'S100'),
    ('read_start_hours_since_expt_start', np.float),
    ('read_start_hours_since_expt_start_qtrhrbucketname', np.float),
    ('read_end_hours_since_expt_start', np.float),
    ('read_end_hours_since_expt_start_qtrhrbucketname', np.float),
    ('read_end_hours_since_read_start', np.float),
    ('read_end_hours_since_read_start_qtrhrbucketname', np.float),
    ('event_count', np.int),
    ('events_per_second', np.float),
    # Analyses
    ('is_analysed', np.int),
    # Analyses/WORKFLOWNAME_000
    ('name', 'S100'),
    ('time_stamp', 'S100'),
    ('nameversion', 'S100'),
    ('has_template', np.int),
    ('has_complement', np.int),
    ('has_2D', np.int),
    # Analyses/WORKFLOWNAME_000/Log
    ('workflow_fullname', 'S100'),
    ('workflow_shortname', 'S100'),
    ('workflow_version', 'S100'),
    # Analyses/WORKFLOWNAME_000/Summary
    ('has_basecall_1d_complement', np.int),
    ('has_basecall_1d_template', np.int),
    ('has_basecall_2d', np.int),
    ('has_hairpin_align', np.int),
    ('has_post_processing_complement', np.int),
    ('has_post_processing_template', np.int),
    ('has_split_hairpin', np.int),
    # Analyses/WORKFLOWNAME_000/Summary/basecall_1d_template
    ('template_num_events', np.int),
    ('template_num_skips', np.int),
    ('template_num_stays', np.int),
    ('template_called_events', np.int),
    ('template_mean_qscore', np.float),
    ('template_strand_score', np.float),
    # Analyses/WORKFLOWNAME_000/Summary/post_processing_template
    ('template_num_merged_events', np.int),
    ('template_raw_events', np.int),
    # Analyses/WORKFLOWNAME_000/Summary/basecall_1d_complement
    ('complement_num_events', np.int),
    ('complement_num_skips', np.int),
    ('complement_num_stays', np.int),
    ('complement_called_events', np.int),
    ('complement_mean_qscore', np.float),
    ('complement_strand_score', np.float),
    # Analyses/WORKFLOWNAME_000/Summary/post_processing_complement
    ('complement_num_merged_events', np.int),
    ('complement_raw_events', np.int),
    # Analyses/WORKFLOWNAME_000/Summary/basecall_2d
    ('twod_mean_qscore', np.float),
    ('twod_sequence_length', np.float),
    # Analyses/WORKFLOWNAME_000/Summary/hairpin_align
    ('hairpin_align_complement_end', np.int),
    ('hairpin_align_complement_start', np.int),
    ('hairpin_align_num_complement', np.int),
    ('hairpin_align_num_template', np.int),
    ('hairpin_align_template_end', np.int),
    ('hairpin_align_template_start', np.int),
    # Analyses/WORKFLOW_000/Summary/split_hairpin
    ('split_hairpin_duration_comp', np.float),
    ('split_hairpin_duration_temp', np.float),
    ('split_hairpin_hairpin_len', np.int),
    ('split_hairpin_num_comp', np.int),
    ('split_hairpin_num_events', np.int),
    ('split_hairpin_num_temp', np.int),
    ('split_hairpin_split_index', np.int),
    # Analyses/WORKFLOWNAME_000/Configuration/general
    ('min_events', np.int),
    ('max_events', np.int),
    ('model_type', 'S50'),
    ('workflow_name', 'S100'),
    ('workflow_script', 'S100'),
    # Analyses/WORKFLOW_NAME/BaseCalled_template
    ('template_events_start_time', np.float),
    ('template_events_duration', np.float),
    ('template_fastq_seqlen', np.int),
    ('template_fastq_bqlen', np.int),
    ('template_fastq_bqmean', np.float),
    ('template_fastq_bqmedian', np.float),
    ('template_numN', np.int),
    # Analyses/WORKFLOW_NAME/BaseCalled_complement
    ('complement_events_start_time', np.float),
    ('complement_events_duration', np.float),
    ('complement_fastq_seqlen', np.int),
    ('complement_fastq_bqlen', np.int),
    ('complement_fastq_bqmean', np.float),
    ('complement_fastq_bqmedian', np.float),
    ('complement_numN', np.int),
    # Analyses/WORKFLOW_NAME/BaseCalled_2D
    # XXXX - Still need to work out how to extract the Alignment from here - template, complement, kmer
    ('twod_fastq_seqlen', np.int),
    ('twod_fastq_bqlen', np.int),
    ('twod_fastq_bqmean', np.float),
    ('twod_fastq_bqmedian', np.float),
    ('twod_numN', np.int),
    ('twod_readclass', 'S16')
]
_call_dtype = None      # The dtype with all the eps_b1-10 values as well
_call_dtype = _rawcall_dtype[:]

# ============================================================================ #
# Program usage                                                                #
# ============================================================================ #

def Initialise():
    'Process command-line arguments, set up output directories and define output filenames.'

    global _progdir, _progname, _args
    _progdir = os.path.dirname(os.path.realpath(sys.argv[0]))
    _progname = os.path.basename(os.path.realpath(sys.argv[0]))
    progexamples = []
    parser = argparse.ArgumentParser(description=_progdesc, \
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--incallsdir', metavar='str', dest='incallsdir', \
        type=str, default=None, help='Input directory containing the raw .fast5 files of event data, or the filtered/the_rest or fail/pass subdirs of fast5 files', required=True)
    parser.add_argument('--inreadstatssimpfile', metavar='sr', dest='inreadstatssimpfile', \
        type=str, default=None, help='', required=True)
    parser.add_argument('--outstatdir', metavar='str', dest='outstatdir', \
        type=str, default=None, help='Output directory for .txt files of read statistics', required=True)
    parser.add_argument('--outprefix', metavar='str', dest='outprefix', \
        type=str, default=None, help='All output files will start with this prefix', required=True)
    parser.add_argument('--overwrite', action='store_true', dest='overwrite', \
        help='Force recomputation of all output files, overwriting existing as necessary', required=False)
    parser.add_argument('--debug', action='store_true', dest='debug', \
        help='Print verbose diagnostics (value must be 0 or 1)', required=False)
    _args = parser.parse_args()

    try:
        if not os.path.exists(os.path.expandvars(_args.outstatdir)):
            os.makedirs(os.path.expandvars(_args.outstatdir))
    except:
        pass
    if not os.path.exists(os.path.expandvars(_args.outstatdir)):
        sys.stderr.write('Erro: Failed to create output directory ({0})'.format(_args.outstatdir))
        sys.exit(_ErrorDirCreate)

    global _single_path
    _single_path = os.path.join(_args.outstatdir, _args.outprefix + '_callstats.txt')

def Read_ReadStatsSimpFile():
    ''
    global _RS
    _RS = dict(np.loadtxt(_args.inreadstatssimpfile, dtype=_readstatssimp_dtype, skiprows=1))

def Get_CallStatsFromFile(fast5_path, twod_readclass):
    'Parse the fast5 analysed file and return a np.array of the statistics.'
    fast5_filename = os.path.basename(fast5_path)
    run_number = '_'.join(fast5_filename.split('_')[-5:-3])
    file_number = fast5_filename.split('_')[-2].replace('file', '')
    if os.stat(fast5_path)[stat.ST_SIZE] == 0:
        sys.stderr.write('Erro: analysed fast5 file has length of zero - ignoring ({0})\n'.format(fast5_path))
        return None
    try:
        hdf = h5py.File(fast5_path, 'r')
    except:
        sys.stderr.write('Erro: Failed to open file with h5py.File for unknown reason - file must be corrupt ({0})\n'.format(fast5_path))
        return None
    try:
        read_number = [int(x.split('_')[1]) for x in hdf['Analyses/EventDetection_000/Reads'].keys() if x.startswith('Read_')][0]
    except:
        sys.stderr.write('Erro: Failed to get read_number from EventDetection_000 section - file must be corrupt ({0})\n'.format(fast5_path))
        return None

    try:
        # Old data
        key = 'UniqueGlobalKey/read_id'
        channel_number = int(hdf[key].attrs['channel_number'])
        read_number2 = int(hdf[key].attrs['read_number'])
        digitisation = 0.0
        offset = 0.0
        range = 0.0
        sampling_rate = _RS[fast5_filename] if fast5_filename in _RS.keys() else _RS[_RS.keys()[0]]
    except:
        # New data
        key = 'UniqueGlobalKey/channel_id'
        channel_number = int(hdf[key].attrs['channel_number'])
        digitisation = float(hdf[key].attrs['digitisation'])
        offset = float(hdf[key].attrs['offset'])
        range = float(hdf[key].attrs['range'])
        sampling_rate = float(hdf[key].attrs['sampling_rate'])

    # Information from the raw events fast5 files needed to convert times.

    key = 'UniqueGlobalKey/tracking_id'
    asic_id = hdf[key].attrs['asic_id']
    asic_temp = float(hdf[key].attrs['asic_temp'])
    device_id = hdf[key].attrs['device_id']
    exp_script_purpose = hdf[key].attrs['exp_script_purpose']
    exp_start_time = int(hdf[key].attrs['exp_start_time'])
    exp_start_localtime_iso = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(exp_start_time))
    flow_cell_id = hdf[key].attrs['flow_cell_id']
    heatsink_temp = float(hdf[key].attrs['heatsink_temp'])
    run_id = hdf[key].attrs['run_id']
    version_name = hdf[key].attrs['version_name']

    key = 'Sequences/Meta'
    numerical_encoding = ''.join(hdf[key].attrs['numerical_encoding'])
    precision = ''.join(hdf[key].attrs['precision'])
    tool = ''.join(hdf[key].attrs['tool'])
    version = ''.join(hdf[key].attrs['version'])

    key = 'Analyses/EventDetection_000/Reads/Read_{0}'.format(read_number)
    keyevents = 'Analyses/EventDetection_000/Reads/Read_{0}/Events'.format(read_number)
    try:
        # Old data
        duration = float(hdf[keyevents].attrs['duration'])
        start_time = float(hdf[keyevents].attrs['start_time'])
    except:
        # New data
        duration = float(hdf[key].attrs['duration']) # Duration is number of 1/sampling_rate-th seconds
        start_time = float(hdf[key].attrs['start_time']) # start_time is number of 1/sampling-rate-th seconds since exp_start_time

    #end_mux = int(hdf[key].attrs['end_mux'])
    ##read_number = hdf[key].attrs['read_number']
    #start_mux = int(hdf[key].attrs['start_mux'])
    try:
        read_start_seconds = exp_start_time + start_time
        read_end_seconds = exp_start_time + (start_time + duration)
        read_start_localtime_iso = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(read_start_seconds))
        read_end_localtime_iso = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(read_end_seconds))
        read_start_hours_since_expt_start = hdf[keyevents][()][0][2] / sampling_rate / 3600.0
        read_start_hours_since_expt_start_qtrhrbucketname = '{0:.2f}'.format(math.ceil(read_start_hours_since_expt_start / 0.25) * 0.25)
        read_end_hours_since_expt_start = hdf[keyevents][()][-1][2] / sampling_rate / 3600.0
        read_end_hours_since_expt_start_qtrhrbucketname = '{0:.2f}'.format(math.ceil(read_end_hours_since_expt_start / 0.25) * 0.25)
        read_end_hours_since_read_start = (hdf[keyevents][()][-1][2] - hdf[keyevents][()][0][2]) / sampling_rate / 3600.0
        read_end_hours_since_read_start_qtrhrbucketname = '{0:.2f}'.format(math.ceil(read_end_hours_since_read_start / 0.25) * 0.25)

        key = 'Analyses/EventDetection_000/Reads/Read_{0}/Events'.format(read_number)
        event_count = len(hdf[keyevents][()])
        events_per_second = (event_count / duration) if (duration and event_count > 100) else 0.0
    except:
        sys.stderr.write('Erro: Failed to extract times from event array, discarding this read and continuing ({0})'.format(fast5_path))
        return None

    is_analysed = 'Analyses' in hdf.keys()
    analyses_workflow = [x for x in hdf['Analyses'].keys() if x.startswith('Basecall')][0]

    key = 'Analyses/{0}'.format(analyses_workflow)
    name = hdf[key].attrs['name']		# metrichor
    time_stamp = hdf[key].attrs['time_stamp']	# YYYY-Mmm-DD HH:MM:SS
    nameversion = hdf[key].attrs['version']	# metrichor version (?)
    has_template = 'BaseCalled_template' in hdf[key].keys()
    has_complement = 'BaseCalled_complement' in hdf[key].keys()
    has_2D = 'BaseCalled_2D' in hdf[key].keys()

    try:
        tmp_log_text = hdf['Analyses/{0}/Log'.format(analyses_workflow)][()].split('\n')
        workflow_fullname = ' '.join(tmp_log_text[0].split(' ')[2:])
        workflow_shortname = tmp_log_text[1].split(' ')[4].strip('.')
        workflow_version = tmp_log_text[0].split('(')[1].split(')')[0].split(' ')[1]
    except:
        workflow_fullname = ''
        workflow_shortname = ''
        workflow_version = ''

    key = 'Analyses/{0}/Summary'.format(analyses_workflow)
    has_basecall_1d_complement = 'basecall_1d_complement' in hdf[key].keys()
    has_basecall_1d_template = 'basecall_1d_template' in hdf[key].keys()
    has_basecall_2d = 'basecall_2d' in hdf[key].keys()
    has_hairpin_align = 'hairpin_align' in hdf[key].keys()
    has_post_processing_complement = 'post_processing_complement' in hdf[key].keys()
    has_post_processing_template = 'post_processing_template' in hdf[key].keys()
    has_split_hairpin = 'split_hairpin' in hdf[key].keys()

    if has_basecall_1d_template:
        key = 'Analyses/{0}/Summary/basecall_1d_template'.format(analyses_workflow)
        template_num_events = hdf[key].attrs['num_events']
        template_num_skips = hdf[key].attrs['num_skips']
        template_num_stays = hdf[key].attrs['num_stays']
        template_called_events = hdf[key].attrs['called_events']
        template_mean_qscore = hdf[key].attrs['mean_qscore']
        template_strand_score = hdf[key].attrs['strand_score']
    else:
        template_num_events = 0
        template_num_skips = 0
        template_num_stays = 0
        template_called_events = 0
        template_mean_qscore = 0.0
        template_strand_score = 0.0
    if has_post_processing_template:
        key = 'Analyses/{0}/Summary/post_processing_template'.format(analyses_workflow)
        template_num_merged_events = hdf[key].attrs['num_merged_events']
        template_raw_events = hdf[key].attrs['num_raw_events']
    else:
        template_num_merged_events = 0
        template_raw_events = 0

    if has_basecall_1d_complement:
        key = 'Analyses/{0}/Summary/basecall_1d_complement'.format(analyses_workflow)
        complement_num_events = hdf[key].attrs['num_events']
        complement_num_skips = hdf[key].attrs['num_skips']
        complement_num_stays = hdf[key].attrs['num_stays']
        complement_called_events = hdf[key].attrs['called_events']
        complement_mean_qscore = hdf[key].attrs['mean_qscore']
        complement_strand_score = hdf[key].attrs['strand_score']
    else:
        complement_num_events = 0
        complement_num_skips = 0
        complement_num_stays = 0
        complement_called_events = 0
        complement_mean_qscore = 0.0
        complement_strand_score = 0.0
    if has_post_processing_complement:
        key = 'Analyses/{0}/Summary/post_processing_complement'.format(analyses_workflow)
        complement_num_merged_events = hdf[key].attrs['num_merged_events']
        complement_raw_events = hdf[key].attrs['num_raw_events']
    else:
        complement_num_merged_events = 0
        complement_raw_events = 0

    if has_basecall_2d:
        key = 'Analyses/{0}/Summary/basecall_2d'.format(analyses_workflow)
        twod_mean_qscore = hdf[key].attrs['mean_qscore']
        twod_sequence_length = hdf[key].attrs['sequence_length']
    else:
        twod_mean_qscore = 0.0
        twod_sequence_length = 0

    if has_hairpin_align:
        key = 'Analyses/{0}/Summary/hairpin_align'.format(analyses_workflow)
        hairpin_align_complement_end = hdf[key].attrs['complement_end']
        hairpin_align_complement_start = hdf[key].attrs['complement_start']
        hairpin_align_num_complement = hdf[key].attrs['num_complement']
        hairpin_align_num_template = hdf[key].attrs['num_template']
        hairpin_align_template_end = hdf[key].attrs['template_end']
        hairpin_align_template_start = hdf[key].attrs['complement_end']
    else:
        hairpin_align_complement_end = 0
        hairpin_align_complement_start = 0
        hairpin_align_num_complement = 0
        hairpin_align_num_template = 0
        hairpin_align_template_end = 0
        hairpin_align_template_start = 0

    if has_split_hairpin:
        key = 'Analyses/{0}/Summary/split_hairpin'.format(analyses_workflow)
        split_hairpin_duration_comp = hdf[key].attrs['duration_comp']
        split_hairpin_duration_temp = hdf[key].attrs['duration_temp']
        try:
            split_hairpin_hairpin_len = hdf[key].attrs['hairpin_len']
        except:
            try:
                split_hairpin_hairpin_len = hdf[key].attrs['hairpin_events']
            except:
                split_hairpin_hairpin_len = 0
        split_hairpin_num_comp = hdf[key].attrs['num_comp']
        split_hairpin_num_events = hdf[key].attrs['num_events']
        split_hairpin_num_temp = hdf[key].attrs['num_temp']
        split_hairpin_split_index = hdf[key].attrs['split_index']
    else:
        split_hairpin_duration_comp = 0.0
        split_hairpin_duration_temp = 0.0
        split_hairpin_hairpin_len = 0
        split_hairpin_num_comp = 0
        split_hairpin_num_events = 0
        split_hairpin_num_temp = 0
        split_hairpin_split_index = 0

    if 'Configuration' in hdf['Analyses/{0}'.format(analyses_workflow)].keys():
        key = 'Analyses/{0}/Configuration/general'.format(analyses_workflow)
        min_events = hdf[key].attrs['min_events']
        max_events = hdf[key].attrs['max_events']
        try:
            model_type = hdf[key].attrs['model_type']
        except:
            model_type = ''
        try:
            workflow_name = hdf[key].attrs['workflow_name']
        except:
            workflow_name =''
        try:
            workflow_script = hdf[key].attrs['workflow_script']
        except:
            workflow_script = ''
    else:
        min_events = 0
        max_events = 0
        model_type = ''
        workflow_name = ''
        workflow_script = ''

    if has_template:
        key = 'Analyses/{0}/BaseCalled_template/Events'.format(analyses_workflow)
        template_events_start_time = hdf[key].attrs['start_time']
        template_events_duration = hdf[key].attrs['duration']
        key = 'Analyses/{0}/BaseCalled_template/Fastq'.format(analyses_workflow)
        template_fastqL = str(hdf['Analyses/{0}/BaseCalled_template/Fastq'.format(analyses_workflow)][()]).split('\n')
        template_fastq_seqlen = len(template_fastqL[1])
        template_fastq_bqlen = len(template_fastqL[3])
        A = np.array([ord(x)-33 for x in str(hdf['Analyses/{0}/BaseCalled_template/Fastq'.format(analyses_workflow)][()]).split('\n')[3]])
        template_fastq_bqmean = np.mean(A)
        template_fastq_bqmedian = np.median(A)
        template_numN = template_fastqL[1].upper().count('N')
    else:
        template_events_start_time = 0.0
        template_events_duration = 0.0
        template_fastq_seqlen = 0
        template_fastq_bqlen = 0
        template_fastq_bqmean = 0.0
        template_fastq_bqmedian = 0.0
        template_numN = 0

    if has_complement:
        key = 'Analyses/{0}/BaseCalled_complement/Events'.format(analyses_workflow)
        complement_events_start_time = hdf[key].attrs['start_time']
        complement_events_duration = hdf[key].attrs['duration']
        key = 'Analyses/{0}/BaseCalled_complement/Fastq'.format(analyses_workflow)
        A = np.array([ord(x)-33 for x in str(hdf['Analyses/{0}/BaseCalled_complement/Fastq'.format(analyses_workflow)][()]).split('\n')[3]])
        complement_fastqL = str(hdf['Analyses/{0}/BaseCalled_complement/Fastq'.format(analyses_workflow)][()]).split('\n')
        complement_fastq_seqlen = len(complement_fastqL[1])
        complement_fastq_bqlen = len(complement_fastqL[3])
        complement_fastq_bqmean = np.mean(A)
        complement_fastq_bqmedian = np.median(A)
        complement_numN = complement_fastqL[1].upper().count('N')
    else:
        complement_events_start_time = 0.0
        complement_events_duration = 0.0
        complement_fastq_seqlen = 0
        complement_fastq_bqlen = 0
        complement_fastq_bqmean = 0.0
        complement_fastq_bqmedian = 0.0
        complement_numN = 0

    if has_2D:
        key = 'Analyses/{0}/BaseCalled_2d/Fastq'.format(analyses_workflow)
        A = np.array([ord(x)-33 for x in str(hdf['Analyses/{0}/BaseCalled_2D/Fastq'.format(analyses_workflow)][()]).split('\n')[3]])
        twod_fastqL = str(hdf['Analyses/{0}/BaseCalled_2D/Fastq'.format(analyses_workflow)][()]).split('\n')
        twod_fastq_seqlen = len(twod_fastqL[1])
        twod_fastq_bqlen = len(twod_fastqL[3])
        twod_fastq_bqmean = np.mean(A)
        twod_fastq_bqmedian = np.median(A)
        twod_numN = twod_fastqL[1].upper().count('N')
    else:
        twod_fastq_seqlen = 0
        twod_fastq_bqlen = 0
        twod_fastq_bqmean = 0.0
        twod_fastq_bqmedian = 0.0
        twod_numN = 0

    call_summary = np.array(
        (
        fast5_filename,
        run_number,
        file_number,
        channel_number,
        read_number,
        sampling_rate,
        asic_id,
        asic_temp,
        device_id,
        exp_script_purpose,
        exp_start_time,
        exp_start_localtime_iso,
        flow_cell_id,
        heatsink_temp,
        run_id,
        version_name,
        numerical_encoding,
        precision,
        tool,
        version,
        duration,
        start_time,
        read_start_localtime_iso,
        read_end_localtime_iso,
        read_start_hours_since_expt_start,
        read_start_hours_since_expt_start_qtrhrbucketname,
        read_end_hours_since_expt_start,
        read_end_hours_since_expt_start_qtrhrbucketname,
        read_end_hours_since_read_start,
        read_end_hours_since_read_start_qtrhrbucketname,
        event_count,
        events_per_second,
        is_analysed,
        name,
        time_stamp,
        nameversion,
        has_template,
        has_complement,
        has_2D,
        workflow_fullname,
        workflow_shortname,
        workflow_version,
        has_basecall_1d_complement,
        has_basecall_1d_template,
        has_basecall_2d,
        has_hairpin_align,
        has_post_processing_complement,
        has_post_processing_template,
        has_split_hairpin,
        template_num_events,
        template_num_skips,
        template_num_stays,
        template_called_events,
        template_mean_qscore,
        template_strand_score,
        template_num_merged_events,
        template_raw_events,
        complement_num_events,
        complement_num_skips,
        complement_num_stays,
        complement_called_events,
        complement_mean_qscore,
        complement_strand_score,
        complement_num_merged_events,
        complement_raw_events,
        twod_mean_qscore,
        twod_sequence_length,
        hairpin_align_complement_end,
        hairpin_align_complement_start,
        hairpin_align_num_complement,
        hairpin_align_num_template,
        hairpin_align_template_end,
        hairpin_align_template_start,
        split_hairpin_duration_comp,
        split_hairpin_duration_temp,
        split_hairpin_hairpin_len,
        split_hairpin_num_comp,
        split_hairpin_num_events,
        split_hairpin_num_temp,
        split_hairpin_split_index,
        min_events,
        max_events,
        model_type,
        workflow_name,
        workflow_script,
        template_events_start_time,
        template_events_duration,
        template_fastq_seqlen,
        template_fastq_bqlen,
        template_fastq_bqmean,
        template_fastq_bqmedian,
        template_numN,
        complement_events_start_time,
        complement_events_duration,
        complement_fastq_seqlen,
        complement_fastq_bqlen,
        complement_fastq_bqmean,
        complement_fastq_bqmedian,
        complement_numN,
        twod_fastq_seqlen,
        twod_fastq_bqlen,
        twod_fastq_bqmean,
        twod_fastq_bqmedian,
        twod_numN,
        twod_readclass
        ),
        dtype=_rawcall_dtype
    )
    hdf.close()
    row = list(np.atleast_1d(call_summary).tolist()[0])
    return row

def CallStats_FileToArray():
    'Read entire callstats.txt file into globals _C and _Cfile.'
    global _C, _Cfile
    _C = np.loadtxt(_single_path, dtype=_call_dtype, delimiter='\t', skiprows=1)
    _Cfile = {}
    for i in range(0, len(_C)):
        _Cfile[_C[i][0]] = 1

def CallStats_Fast5ToFile(indir, filelist, readclass):
    'Write all existing data from _C to callstats.txt, then append rows from new files in filelist.'

    # If callstats.txt file already exists, rename it so it is can be accessed later if necessary.'
    if os.path.exists(_single_path):
        new_single_path = _single_path.replace('.txt', '_' + time.strftime('%Y%m%d-%H%M%S', time.localtime()) + '.txt')
        try:
            os.rename(_single_path, new_single_path)
        except:
            sys.stderr.write('Error: Failed to make backup existing callstats.txt file ({0} to {1})\n'.format(_single_path, new_single_path))
            sys.exit(_ErrorFileRename)

    # Open the output readstats.txt file.
    lenfilelist = len(filelist)
    single_fp = open(_single_path, 'w')
    single_fp.write('{0}\n'.format('\t'.join([x[0] for x in _call_dtype])))
    # dump all the existing data from _C into it.
    if _C is not None:
        for i in range(0, len(_C)):
            row = [str(x) for x in _C[i]]
            single_fp.write('{0}\n'.format('\t'.join(row)))
        single_fp.flush()
    for j, file in enumerate(filelist):
        if _args.debug and (j+1) >= 500:
            break
        thing_path = os.path.join(indir, file)
        if _Cfile.has_key(file):
            continue
        row = Get_CallStatsFromFile(thing_path, readclass)
        if row is not None:
            single_fp.write('{0}\n'.format('\t'.join([str(x) for x in row])))
            if (((j+1) % 100) == 0) or ((j+1) == lenfilelist):
                sys.stdout.write('{0} Info: Processing {1}-th of {2} additional fast5 file\n'.format(
                    time.strftime('%Y%m%d-%H%M%S', time.localtime()), j+1, lenfilelist))
                single_fp.flush()
        else:
            sys.stdout.write('{0} Info: Processing failed for {1}-th additional fast5 file ({2})\n'.format(
                time.strftime('%Y%m%d-%H%M%S', time.localtime()), j+1, thing_path))
    # And close the callstats.txt file.
    if single_fp:
        single_fp.close()

def Extract_CallStats_FromDir(incallsdir, readclass):
    'Iterate through all reads/downloads/prefix.fast5 files and output a single callstats/callstats.txt file.'

    # Set up IO directories.
    outstatdir = os.path.expandvars(_args.outstatdir)

    # Initialise global variables.
    global _C, _Cfile

    # Get current list of raw fast5 files to be processed.
    files = [x for x in os.listdir(incallsdir) if os.path.isfile(os.path.join(incallsdir, x)) and x.endswith('.fast5')]
    files.sort()

    # If callstats.txt file exists ...
    #     check if callstats.txt file is complete and go to next section if it is.
    # ... if callstats.txt does not exist, update with info from all basecalled fast5 files.
    if os.path.exists(_single_path):
        CallStats_FileToArray()
        numrecords_before = len(_C)
        new_set = set(files)
        cur_set = set(_Cfile.keys())
        if new_set.issubset(cur_set) and not _args.overwrite:
            sys.stdout.write('Info: callstats.txt already contains data from calldir, no need to compute ({0})\n'.format(incallsdir))
            numrecords_after = len(_C)
        else:
            if not _args.overwrite:
                sys.stdout.write('Info: callstats.txt exists but may not be complete, starting update now ({0})\n'.format(_single_path))
            else:
                sys.stdout.write('Info: callstats.txt exists but overwrite flag is forcing recomputation\n')
                _C = None
                _Cfile = {}
            CallStats_Fast5ToFile(incallsdir, files, readclass)
            CallStats_FileToArray()
            numrecords_after = len(_C)
    else:
        numrecords_before = 0
        sys.stdout.write('Info: callstats.txt does not exist, starting iteration over {0} fast5 files now ({1})\n'.format(len(files), _single_path))
        CallStats_Fast5ToFile(incallsdir, files, readclass)
        CallStats_FileToArray()
        numrecords_after = len(_C)
    global _callstatsfile_was_updated
    _callstatsfile_was_updated = (numrecords_before < numrecords_after) or _args.overwrite
    sys.stdout.write('{0} Info: callstats.txt updated with {1} new records\n'.format(time.strftime('%Y%m%d-%H%M%S', time.localtime()), numrecords_after - numrecords_before))

    return

    for thing in os.listdir(incallsdir):
        thing_path = os.path.join(incallsdir, thing)
        if not os.path.isfile(thing_path) or not thing.endswith('.fast5'):
            continue
        thing_stem = thing[:-6]
        stat_path = os.path.join(_args.outstatdir, thing_stem + '_callstats.txt')
        if os.path.exists(stat_path):
            continue
        row = Get_CallStatsFromFile(thing_path, readclass)
        if row is not None:
            with open(stat_path, 'w') as out_fp:
                out_fp.write('{0}\n'.format('\t'.join([x[0] for x in _call_dtype])))
                out_fp.write('{0}\n'.format('\t'.join([str(x) for x in row])))

def Extract_CallStats():
    '''
    Iterate through all reads/downloads/subdir/prefix.fast5 files where subdir could be
    filtered, the_rest, fail or pass, and output a single callstats/callstats.txt file.
    '''
    global _C, _Cfile
    _C = None
    _Cfile = {}

    if os.path.exists(os.path.join(_args.incallsdir, 'filtered')):
        Extract_CallStats_FromDir(os.path.join(_args.incallsdir, 'filtered'), 'pass')
        Extract_CallStats_FromDir(os.path.join(_args.incallsdir, 'the_rest'), 'fail')
    elif os.path.exists(os.path.join(_args.incallsdir, 'pass')):
        Extract_CallStats_FromDir(os.path.join(_args.incallsdir, 'pass'), 'pass')
        Extract_CallStats_FromDir(os.path.join(_args.incallsdir, 'fail'), 'fail')
    else:
        Extract_CallStats_FromDir(os.path.join(_args.incallsdir))

# ============================================================================ #
# Main                                                                         #
# ============================================================================ #

if __name__ == '__main__':
    Initialise()
    Read_ReadStatsSimpFile()
    Extract_CallStats()

# ============================================================================ #
