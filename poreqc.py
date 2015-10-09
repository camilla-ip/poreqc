#!/usr/bin/env python

# ============================================================================ #
# poreqc.py                                                                    #
# ============================================================================ #
"""
PoreQC MinION Sequencing Report

Produces a PDF and HTML report based on either the intermediate or final output
data from an Oxford Nanopore MinION device.

Input:
- a file of meta-data related to sample provenance, library preparation and
  suitable reference for the target genome(s) with which to do the QC.
- directory containing a copy of the basecalled and/or pre-basecalled .fast5 files
  produced by the run so far.

Output:
- a directory containing one text file of tabular data required for each plot
- a directory of plots with suitable size and resolution for inclusion in a PDF
- a directory of similar plots suitable for an HTML report to be displayed on a laptop [NOT IMPL]
- a printable PDF version of the report, generated via LaTeX
- an HTML version of the report that refers to images in a sub-directory [NOT IMPL]

Notes:
1. The program was designed to be run following an rsync of the sequencing data
   from the sequencing laptop. A special 'data status file' will indicate whether
   the run has completed or is still ongoing.
2. It is expected that the output directory for the HTML version of the report
   will be a directory which is accessed by a web server so the laboratory people
   running the MinION can view the current status of the data in real-time.
"""
# ============================================================================ #
# Camilla Ip                                                                   #
# camilla.ip@well.ox.ac.uk                                                     #
# July 2014                                                                    #
# ============================================================================ #

# ============================================================================ #
# Import Modules                                                               #
# ============================================================================ #

import multiprocessing as mp
import subprocess as sp
import argparse, datetime, h5py, math, numpy as np, os, random, re, shlex, stat, sys, time
from Bio import SeqIO
from Bio.SeqUtils import GC
from StringIO import StringIO
from subprocess import PIPE

sys.path.insert(0, '/'.join(os.path.dirname(os.path.realpath(sys.argv[0])).split('/')[:-2] + ['bin']))
sys.path.insert(0, os.path.dirname(os.path.realpath(sys.argv[0])))

# ============================================================================ #
# Global variables                                                             #
# ============================================================================ #

_progdir = None
_progname = None
_intdir = None
_version = '0.2.10'
_progdesc = 'poreqc sequencing run report v{0}'.format(_version)

_ErrorFileMissing = 23

_currentdatetime = None		# The current system date/time.

_S = None			# Current status of the sequencing run.
				# Must have a single line containing the word "Running" or "Finished".
_metadatafile_fields = [
    'runstartdate',
    'hostname',
    'flowcellid',
    'samplenameshort',
    'outfileprefix',
    'experimentname',
    'experimentaims',
    'experimentdesc',
    'samplenamefull',
    'species',
    'chemistryversion',
    'programlength_hrs',
    'minknowversion',
    'maprefid',
    'mindepthrequired'
]
_M = {}         # Metadata dictionary _M[field] = value
sectionD = {}	# Pre-computed values for the tabular parts of the report.

_R = None	# The SeqIO object holding the reference genome. If no reference
		# is provided, this value will continue to have the value of None.

_outpath = {}	# Set up all the output paths at the beginning
_outfp = {}	# File pointers to the outpaths

_ML = []	# Ordered list of metadata sections and fields _ML[section] = [field, ...]

_read_dtype = [
    # From filename _chN_fileN.fast5
    ('fast5_filename', 'S100'),
    ('run_number', 'S100'),
    ('file_number', np.int),
    # /Key/read_id
    ('channel_number', np.int),
    ('read_number', np.int),
    # /Key/tracking_id
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
    # /Reads/Read_READ_NUMBER
    ('event_count', np.int),
    ('median_before', np.float),
    ('sample_rate', np.float),
    ('start', np.int),
    ('start_mux', np.int),
    # /Sequences/Meta
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
    ('eps_2ndhalf', np.float),
    ('eps_b1', np.float),
    ('eps_b2', np.float),
    ('eps_b3', np.float),
    ('eps_b4', np.float),
    ('eps_b5', np.float),
    ('eps_b6', np.float),
    ('eps_b7', np.float),
    ('eps_b8', np.float),
    ('eps_b9', np.float),
    ('eps_b10', np.float)
]
_E = None	# Event data as a 2D numpy array.

_metrichor_version = 'UNKNOWN'
_workflow_name = 'UNKNOWN'

# This has to be generic enough to be applicable to the DNA (2D) and cDNA (1D) workflows.
_call_dtype = [
    # From filename _chN_fileN.fast5
    ('fast5_filename', 'S100'),
    ('run_number', 'S100'),
    ('file_number', np.int),
    ('read_number', np.int),	# From reads section, but putting it here to be consistent with read_dtype.
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

_bcall_dtype= [
    ('fast5_filename', 'a200'),
    ('ch', np.int),
    ('file', np.int),
    ('read_type_from_log', 'a100'),
    ('channel_number', np.int),
    ('read_number', np.int),
    ('device_id', 'S100'),
    ('asic_id', 'S100'),
    ('flow_cell_id', 'S100'),
    ('asic_temp', np.float),
    ('heatsink_temp', np.float),
    ('exp_start_time', np.int),
    ('exp_start_time_nice', 'S50'),
    ('version_name', 'S100'),
    ('run_id', 'S100'),
    ('is_analysed', np.bool),
    ('event_count', np.int),
    ('is_base_called_2d', np.bool),
    ('has_template_read', np.bool),
    ('has_complement_read', np.bool),
    ('has_2D_read', np.bool),
    ('read_length_template', np.int),
    ('read_length_complement', np.int),
    ('read_length_2D', np.int),
    ('start_time_read', np.float),
    ('start_time_read_nice', 'S50'),
    ('duration_read', np.float),
    ('end_time_read_nice', 'S50'),
    ('start_time_template', np.float),
    ('start_time_template_nice', 'S50'),
    ('duration_template', np.float),
    ('end_time_template_nice', 'S50'),
    ('start_time_complement', np.float),
    ('start_time_complement_nice', 'S50'),
    ('duration_complement', np.float),
    ('end_time_complement_nice', 'S50'),
    ('basecall_called_events_template', np.int),
    ('basecall_num_events_template', np.int),
    ('basecall_num_skips_template', np.int),
    ('basecall_num_stays_template', np.int),
    ('basecall_mean_qscore_template', np.float),
    ('basecall_strand_score_template', np.float),
    ('basecall_called_events_complement', np.int),
    ('basecall_num_events_complement', np.int),
    ('basecall_num_skips_complement', np.int),
    ('basecall_num_stays_complement', np.int),
    ('basecall_mean_qscore_complement', np.float),
    ('basecall_strand_score_complement', np.float),
    ('basecall_mean_qscore_2D', np.float),
    ('basecall_sequence_length_2D', np.int),
    ('basecall_num_raw_events_template', np.int),
    ('basecall_num_merged_events_template', np.int),
    ('basecall_num_raw_events_complement', np.int),
    ('basecall_num_merged_events_complement', np.int)
]
_B = None	# Bcall data as a 2D numpy array.
_Bread = {}	# _Bread[(channel_number, read_number)] = 1 if that channel/read combination has already been read.

_readlenaccum_dtype = [
    ('readlenmx', np.int),
    ('eventtotal_mbases', np.float)
]

_readstats_dtype = [
    ('ch', np.int),
    ('file', np.int),
    ('channel_number', np.int),
    ('read_number', np.int),
    ('event_count', np.int),
    ('read_length_template', np.int),
    ('read_length_complement', np.int),
    ('read_length_2D', np.int),
    ('duration_read', np.float),
    ('duration_template', np.float),
    ('duration_complement', np.float),
    ('read_start_time', np.float),
    ('read_end_time', np.float),
    ('template_start_time', np.float),
    ('template_end_time', np.float),
    ('complement_start_time', np.float),
    ('complement_end_time', np.float)
]

_channelstats_dtype = [
    ('channel_number', np.int),
    ('num_reads', np.int),
    ('num_template_reads', np.int),
    ('num_complement_reads', np.int),
    ('num_2D_reads', np.int),
    ('sum_event_count_mbases', np.float),
    ('sum_template_lengths_mbases', np.float),
    ('sum_complement_lengths_mbases', np.float),
    ('sum_2D_lengths_mbases', np.float)
]

_timestats_start_time = 0	# Experiment start time in seconds since epoch.
				# This is the baseline time from which all off-set
				# times in _timestats_dtype are computed.
_timestats_dtype = [
    ('ch', np.int),
    ('file', np.int),
    ('channel_number', np.int),
    ('read_number', np.int),
    ('time_bucket_in_hours', np.float),
]

_D = {}		# A dictionary containing all the information retrieved from the .fast5 files

# ============================================================================ #
# Program usage                                                                #
# ============================================================================ #

def Initialise():
    '''
    Read in the command-line arguments. Check for errors. Set up the output file
    paths in _outpath dictionary.
    '''

    # Process command-line arguments.
    global _progdir, _progname, _args
    _progdir = os.path.dirname(os.path.realpath(sys.argv[0]))
    _progname = os.path.basename(os.path.realpath(sys.argv[0]))
    progexamples = [
    ]
    parser = argparse.ArgumentParser(description=_progdesc, \
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--inrundir', metavar='str', dest='inrundir', \
        type=str, default=None, help='Top-level directory containing a copy of the data from a sequencing run', required=True)
    parser.add_argument('--metadatafile', metavar='str', dest='metadatafile', \
        type=str, default='metadata.txt', help='Filename of the metadata file that should be in --outdir', required=False)
    parser.add_argument('--statusfile', metavar='str', dest='statusfile', \
        type=str, default='status.txt', help='Filename of the status file that should be in --outdir (Must contain either "Running" or "Finished")', required=False)
    parser.add_argument('--refdir', metavar='str', dest='refdir', \
        type=str, default='/bsgdata/microbial/data/references', help='Top-level directory for all references ' \
        '(e.g., for wtchgR00000021 you would specify /bsgdata/microbial/data/references to find file ' \
        'wtchgR00000021/wtchgR00000021.fasta)', required=True)
    parser.add_argument('--outdir', metavar='str', dest='outdir', \
        type=str, default=None, help='Write all output to this directory', required=True)
    parser.add_argument('--outprefix', metavar='str', dest='outprefix', \
        type=str, default=None, help='All output files will start with this prefix', required=True)
    parser.add_argument('--bconly', action='store_true', dest='bconly', \
        help='Extract all information from the basecalled fast5 files only', required=False)
    parser.add_argument('--forcereadstep', action='store_true', dest='forcereadstep', \
        help='Force recompute of readstats.txt and associated summary tables', required=False)
    parser.add_argument('--forcecallstep', action='store_true', dest='forcecallstep', \
        help='Force recompute of callstats.txt', required=False)
    parser.add_argument('--forcereportstep', action='store_true', dest='forcereportstep', \
        help='Force recompute of PDF report and associated images and summary statistics', required=False)
    parser.add_argument('--overwrite', action='store_true', dest='overwrite', \
        help='Force recomputation of all output files, overwriting existing as necessary', required=False)
    parser.add_argument('--readlenbucketsz', metavar='int', dest='readlenbucketsz', \
        type=int, default=1000, help='Size of read length buckets on x-axis of MinKNOW-like read yield distribution plot', required=False)
    parser.add_argument('--imgszconst', metavar='int', dest='imgszconst', \
        type=int, default=100, help='Size parameter for images', required=False)
    parser.add_argument('--debug', action='store_true', dest='debug', \
        help='Print verbose diagnostics (value must be 0 or 1)', required=False)

    _args = parser.parse_args()

    # Save the current time which is useful for the report date and output sub-directories later
    global _currentdatetime
    _currentdatetime = time.localtime()

    # Set up the output file paths.
    global _outpath
    outstem = os.path.join(os.path.expandvars(_args.outdir), _args.outprefix)
    #_outpath['subdir'] = os.path.join(os.path.expandvars(_args.outdir), time.strftime('%Y%m%d-%H%M%S', _currentdatetime))
    _outpath['subdir'] = os.path.expandvars(_args.outdir)
    _outpath['readstatdir'] = os.path.expandvars(_args.outdir)
    _outpath['readstatsfile'] = os.path.join(_outpath['readstatdir'], _args.outprefix + '_readstats.txt')
    _outpath['simpfile'] = os.path.join(_outpath['readstatdir'], _args.outprefix + '_readstats_simplified.txt')
    _outpath['callstatdir'] = os.path.expandvars(_args.outdir)
    _outpath['callstatsfile'] = os.path.join(_outpath['callstatdir'], _args.outprefix + '_callstats.txt')
    _outpath['activechannelsfile'] = os.path.join(_outpath['readstatdir'], _args.outprefix + '_activechannels.txt')
    _outpath['filldelayfile'] = os.path.join(_outpath['readstatdir'], _args.outprefix + '_filldelay.txt')
    _outpath['pngpdf'] = os.path.join(_outpath['readstatdir'], 'pngpdf')
    _outpath['tabletxt'] = os.path.join(_outpath['subdir'], _args.outprefix + '_table.txt')
    _outpath['tablecsv'] = os.path.join(_outpath['subdir'], _args.outprefix + '_table.csv')
    _outpath['report'] = os.path.join(_outpath['subdir'], _args.outprefix + '_report_v{0}.tex'.format(_version))
    _outpath['reportpdf'] = os.path.join(_outpath['subdir'], _args.outprefix + '_report_v{0}.pdf'.format(_version))
    _outpath['runparameters'] = os.path.join(_outpath['subdir'], _args.outprefix + 'runparams.txt')
    _outpath['tex2pdflog'] = os.path.join(_outpath['subdir'], _args.outprefix + '_report.log')
    _outpath['tex2pdfaux'] = os.path.join(_outpath['subdir'], _args.outprefix + '_report.aux')

    # Create the output sub-directory
    if not os.path.exists(_outpath['subdir']):
        os.makedirs(_outpath['subdir'])
        if not os.path.exists(_outpath['subdir']):
            sys.stderr.write('Erro: Failed to create output dir ({0})\n'.format(_outpath['subdir']))
            sys.exit(2)
    if not os.path.exists(_outpath['readstatdir']):
        os.makedirs(_outpath['readstatdir'])
        if not os.path.exists(_outpath['readstatdir']):
            sys.stderr.write('Erro: Failed to create output dir ({0})\n'.format(_outpath['readstatdir']))
            sys.exit(2)

    # If everything up-to-date then nothing to do.
    if os.path.exists(_outpath['readstatsfile']) and os.stat(_outpath['readstatsfile']).st_size > 0 and \
       os.path.exists(_outpath['simpfile']) and os.stat(_outpath['simpfile']).st_size > 0 and \
       os.path.exists(_outpath['callstatsfile']) and os.stat(_outpath['callstatsfile']).st_size > 0 and \
       os.path.exists(_outpath['activechannelsfile']) and os.stat(_outpath['activechannelsfile']).st_size > 0 and \
       os.path.exists(_outpath['filldelayfile']) and os.stat(_outpath['filldelayfile']).st_size > 0 and \
       os.path.exists(_outpath['tabletxt']) and os.stat(_outpath['tabletxt']).st_size > 0 and \
       os.path.exists(_outpath['tablecsv']) and os.stat(_outpath['tablecsv']).st_size > 0 and \
       os.path.exists(_outpath['reportpdf']) and os.stat(_outpath['reportpdf']).st_size > 0 and \
    not _args.forcereadstep and \
    not _args.forcecallstep and \
    not _args.forcereportstep:
        sys.stdout.write('Info: Nothing to do - exiting now ({0}, {1})\n'.format(_outpath['tablecsv'], _outpath['pngpdf']))
        sys.exit(0)

# ============================================================================ #
# Generic utility functions                                                    #
# ============================================================================ #

def sys_exec(cmd):
    """
    Execute a command using the subprocess module to trap the
    return value, the stdout string and stderr string.
    """
    proc_handle = sp.Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
    proc_stdout, proc_stderr = proc_handle.communicate()
    proc_returncode = proc_handle.returncode
    return [proc_returncode, proc_stdout, proc_stderr]

# ============================================================================ #
# Read external files                                                          #
# ============================================================================ #

def Read_Status():
    'Read the status code in the statusfile and store it in _S.'
    global _S
    status_path = os.path.join(os.path.expandvars(_args.outdir), _args.statusfile)
    with open(status_path, 'r') as status_fp:
        for line in status_fp:
            line = line.strip()
            if not len(line) or line.startswith('#'):
                continue
            if not _S:
                _S = line

def Read_Metadata():
    '''
    Read the contents of the metadata file into the global _M dictionary.
    If any fields are missing, print an error and exit.
    '''
    global _M
    meta_path = os.path.join(os.path.expandvars(_args.outdir), _args.metadatafile)
    if not os.path.exists(meta_path):
        sys.stderr.write('Erro: metadatafile does not exist ({0})\n'.format(meta_path))
        sys.exit(2)
    with open(meta_path, 'r') as meta_fp:
        for line in meta_fp:
            line = line.strip()
            if not len(line) or line.startswith('#'):
                continue
            L = line.split(',')
            if len(L) < 2:
                sys.stderr.write('Erro: metadata line does not contain two comma-separated fields ({0})\n'.format(line))
                sys.exit(2)
            field = L[0]
            val = ','.join(L[1:])
            if field not in _metadatafile_fields:
                sys.stderr.write('Erro: metadata line contains unrecognised field name ({0})\n'.format(line))
                sys.exit(2)
            _M[field] = val
    fieldL = _M.keys()
    for field in _metadatafile_fields:
        if field not in fieldL:
            sys.stderr.write('Erro: metadata for field is missing in file ({0})\n'.format(field))
            sys.exit(2)
    if _args.debug:
        sys.stdout.write('Info: metadata_file ({0})\n'.format(meta_path))
        for field in _metadatafile_fields:
             sys.stdout.write('Info:     {0},{1}\n'.format(field, _M[field]))

def Read_Reference():
    'Read reference, return the SeqIO object which may contain more than one contig.'
    global _R
    _R = {}
    if _M['maprefid'].lower() == 'none':
        sys.stdout.write('Warn: No reference id specified ({0})\n'.format(_M['maprefid']))
        return
    in_path = os.path.expandvars(os.path.join(_args.refdir, _M['maprefid'], _M['maprefid']+'.fasta'))
    if not os.path.exists(in_path):
        sys.stderr.write('Erro: reference FASTA for mapping does not exist ({0})\n'.format(in_path))
        _R = None
        sys.exit(_ErrorFileMissing)
    in_fp = open(in_path, 'r') if not in_path.endswith('gz') else gzip.open(in_path, 'rb')
    for record in SeqIO.parse(in_fp, 'fasta'):
        _R[record.id] = record
    in_fp.close()
    if _args.debug:
        sys.stdout.write('Info: reference has {0} contig(s) with total length {1} ({2}) \n'.format(
            len(_R.keys()),
            sum([len(_R[id].seq) for id in _R.keys()]),
            in_path))

# ============================================================================ #
# Extract data                                                                 #
# ============================================================================ #

def Extract_NewEventData_OneDir(inreadsdir, readsarecalled):
    '''
    Run extract_readstats.py to get one _readstats.txt file for each /reads/downloads/pass/*.fast5
    and /reads/downloads/fail/*.fast5 file and generate a single readstats.txt file containing all the information.
    '''

    prog = os.path.join(os.path.dirname(os.path.expandvars(sys.argv[0])), 'poreqc_extractreadstats.py')
    cmd = '{prog} --inreadsdir {indir} --outstatdir {outdir} --outprefix {outprefix}{overwrite}{debug}{readsarecalled}'.format(
        prog=prog, indir=inreadsdir, outdir=_outpath['readstatdir'],
        outprefix=_args.outprefix,
        overwrite=' --overwrite' if _args.forcereadstep else '',
        debug=' --debug' if _args.debug else '',
        readsarecalled=' --inreadsarecalled' if readsarecalled else '')
    print 'Info: Running: {0}'.format(cmd)
    retval, retout, reterr = sys_exec(cmd)
    if len(retout):
        sys.stdout.write('{0}{1}'.format(retout, '\n' if retout[-1]!='\n' else ''))
    if len(reterr):
        sys.stdout.write('{0}{1}'.format(reterr, '\n' if reterr[-1]!='\n' else ''))
    sys.stdout.write('Command {0}: returncode {1}\n'.format('succeeded' if retval==0 else 'failed', retval))
    if retval != 0:
        sys.exit(8)

def Extract_NewEventData():
    'If bconly is True, use the pass and fail reads, otherwise use the pre-basecalled reads.'

    if _args.bconly:
        passdir = os.path.join(os.path.expandvars(_args.inrundir), 'reads', 'downloads', 'pass')
        faildir = os.path.join(os.path.expandvars(_args.inrundir), 'reads', 'downloads', 'fail')
        Extract_NewEventData_OneDir(passdir, 1)
        Extract_NewEventData_OneDir(faildir, 1)
    else:
        readsdir = os.path.join(os.path.expandvars(_args.inrundir), 'reads')
        Extract_NewEventData_OneDir(readsdir, 0)

def Read_AllEventData():
    'Read contents of readstats.txt into global array _E.'

    global _E
    _E = np.loadtxt(_outpath['readstatsfile'], dtype=_read_dtype, delimiter='\t', skiprows=1)
    pass

def Extract_NewCalldata():
    '''
    Run extract_callstats.py to get one _callstats.txt file for each /reads/downloads/*.fast5 file
    and generate a single callstats.txt file containing all the information.
    '''

    prog = os.path.join(os.path.dirname(os.path.expandvars(sys.argv[0])), 'poreqc_extractcallstats.py')
    incallsdir = os.path.join(os.path.expandvars(_args.inrundir), 'reads', 'downloads')
    cmd = '{prog} --incallsdir {indir} --inreadstatssimpfile {simpfile} --outstatdir {outdir} --outprefix {outprefix}{overwrite}{debug}'.format(
        prog=prog, indir=incallsdir, simpfile=_outpath['simpfile'], outdir=_outpath['callstatdir'],
        outprefix=_args.outprefix,
        overwrite=' --overwrite' if  _args.forcecallstep else '',
        debug=' --debug' if _args.debug else '')
    print 'Info: Running: {0}'.format(cmd)
    retval, retout, reterr = sys_exec(cmd)
    if len(retout):
        sys.stdout.write('{0}{1}'.format(retout, '\n' if retout[-1]!='\n' else ''))
    if len(reterr):
        sys.stdout.write('{0}{1}'.format(reterr, '\n' if reterr[-1]!='\n' else ''))
    sys.stdout.write('Command {0}: returncode {1}\n'.format('succeeded' if retval==0 else 'failed', retval))
    if retval != 0:
        sys.exit(8)

def Read_AllCallData():
    'Read contents of callstats.txt into global array _C.'

    global _C
    _C = np.loadtxt(_outpath['callstatsfile'], dtype=_call_dtype, delimiter='\t', skiprows=1)
    pass

def FixBq(seq, bq):
    ''
    if len(seq) > len(bq):
        newbq = bq + bq[-1]*(len(seq)-len(bq))
        fixtype = 'extended'
        sys.stdout.write('Warn: {0} read of fast5 created with bq string extended with last quality ({1} to {2} with bq num {3} = char {4})\n'.format(readtype, len(bq), len(newbq), ord(bq[-1])-33, bq[-1]))
    elif len(seq) < len(bq):
        newbq = bq[0:len(seq)]
        sys.stdout.write('Warn: {0} read of fast5 created with truncated bq string ({1} to {2} chars)\n'.format(readtype, len(bq), len(newbq)))
        fixtype = 'truncated'
    else:
        newbq = bq
        fixtype = 'none'
    return newbq, fixtype

# ============================================================================ #
# Create plots                                                                 # 
# ============================================================================ #

def Create_Plots():
    '''
    Create all the plots, one after another, in two sizes: one for embedding
    in the PDF report, and another that can be imported in an HTML page.
    '''

    if os.path.exists(_outpath['reportpdf']) and not  _args.overwrite and not _args.forcereportstep:
        sys.stdout.write('Info: PDF report already up to date, no need to re-generate PNG images\n')
        return
    sys.stdout.write('Info: Generating PNG images for PDF report\n')

    Rprog = 'Rscript'
    plotprog = os.path.join(os.path.dirname(os.path.expandvars(sys.argv[0])), 'poreqc_makeplots.R')
    cmd = '{Rprog} {plotprog} {readstatspath} {callstatspath} {activechannelspath} {filldelaypath} {rundesc} {outdir} {outprefix} {imgszconst}'.format(
        Rprog=Rprog,
        plotprog=plotprog,
        readstatspath=_outpath['readstatsfile'],
        callstatspath=_outpath['callstatsfile'],
        activechannelspath=_outpath['activechannelsfile'],
        filldelaypath=_outpath['filldelayfile'],
        rundesc=_M['outfileprefix'],
        outdir=os.path.expandvars(_args.outdir),
        outprefix=_M['outfileprefix'],
        imgszconst=_args.imgszconst
    )
    sys.stdout.write('Info: Running: {0}\n'.format(cmd))
    retval, retout, reterr = sys_exec(cmd)
    if retval != 0:
        sys.stderr.write('Command failed: {0}\nstdout msg: {1}\nstderr msg: {2}'.format(cmd, retout, reterr))
        sys.exit(8)


# ============================================================================ #
# Create report                                                                #
# ============================================================================ #

def escaped(text):
    new_text = re.sub('%', '\%', text)
    new_text = re.sub('\$', '\\$', new_text)
    new_text = re.sub('_', '\_', new_text)
    return new_text

def Report_Init():
    "Headers for the latex document."

    global _outfp

    _outfp['report'] = open(_outpath['report'], 'w')
    if not _outfp:
        sys.stderr.write('Erro: Failed to open output report path ({0})'.format(_outpath['report']))
        sys.exit(2)

    _outfp['report'].write("\documentclass[a4paper,10pt]{article}\n")
    _outfp['report'].write("\\usepackage{graphicx}\n")
    _outfp['report'].write("\\usepackage[margin=0.2in]{geometry}\n")
    _outfp['report'].write("\\usepackage{times}\n")
    _outfp['report'].write("\\usepackage{subfig}\n")
    _outfp['report'].write("\\usepackage{float}\n")
    _outfp['report'].write("\\usepackage[scaled=0.95]{helvet}\n")
    _outfp['report'].write("\\usepackage{fullpage}\n")
    _outfp['report'].write("\\usepackage{bm}\n")
    _outfp['report'].write("\\usepackage{array}\n")

    _outfp['report'].write("\\usepackage{array}\n")
    _outfp['report'].write("\\newcolumntype{L}[1]{>{\\raggedright\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}}\n")
    _outfp['report'].write("\\newcolumntype{C}[1]{>{\\centering\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}}\n")
    _outfp['report'].write("\\newcolumntype{R}[1]{>{\\raggedleft\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}}\n")

    _outfp['report'].write("\\usepackage[font=small,skip=0pt]{caption}\n")

    _outfp['report'].write("\\usepackage[T1]{fontenc}\n")

    _outfp['report'].write("\\addtolength{\\oddsidemargin}{-.5in}\n")
    _outfp['report'].write("\\addtolength{\\evensidemargin}{-.5in}\n")
    _outfp['report'].write("\\addtolength{\\textwidth}{1.0in}\n")
    _outfp['report'].write("\\addtolength{\\topmargin}{-0.5in}\n")
    _outfp['report'].write("\\addtolength{\\textheight}{1.0in}\n")

    _outfp['report'].write("\\begin{document}\n")
    _outfp['report'].write("\\author{\n")
    _outfp['report'].write("    {0}\\\\\n".format(_progdesc))
    _outfp['report'].write("}\n")
    _outfp['report'].write("\\title{{{0}}}\n".format(escaped(_M['experimentname'])))
    _outfp['report'].write("\\maketitle\n")
    _outfp['report'].write("\n")

def Report_Table(data, caption, small, outfp):
    'Turn an array of arrays, or list of lists, into a nicely-formatted latex table.'
    outfp.write("\\begin{quote}\n")
    outfp.write("\\begin{table}[!htpb]\n")
    #reduce the size of tables which would otherwise be too wide for the document
    if small == True:
        outfp.write("\\small\n")
        outfp.write("\\tabcolsep=0.2cm\n")
    outfp.write("\\centering\n")
    outfp.write("\\begin{tabular}{ R{5cm} L{10cm} }\n")
    for line in data:
        outfp.write("{0}\\\\\n".format("&".join([escaped(str(x)) for x in line])))
    outfp.write("\\end{tabular}\n")
    outfp.write("\\end{table}\n")
    outfp.write("\\end{quote}\n")

def Report_BulletList(lineL, outfp):
    ''
    if not len(lineL):
        return
    outfp.write('\\begin{itemize}\n')
    for line in lineL:
        outfp.write('  \\item {0}\n'.format(line))
    outfp.write('\\end{itemize}\n')

def Report_Image(image1, captiontext, captiondesc, outfp):
    'Output one image, centred.'
    outfp.write("\\begin{figure}[H]\n")
    outfp.write("\\captionsetup{justification=raggedright, singlelinecheck=false}\n")
    outfp.write("    \\includegraphics[width=\\textwidth]{{{image}}}\n".format(image=image1))
    #outfp.write("    \\caption{{{caption}}}{{{desc}}}\n".format(caption=captiontext, desc=captiondesc))
    outfp.write("    \\caption{{{caption}}}{{{desc}}}\n".format(caption=captiontext, desc=''))
    outfp.write("\\end{figure}\n")

def Report_Images_Sidebyside(image1, subcaption1, image2, subcaption2, maincaption, outfp):
    'Print 2 images side-by-side to the latex document'
    outfp.write("\\begin{quote}\n")
    outfp.write("\\begin{figure}[H]\n")
    outfp.write("\\subfloat[{caption}]{{\includegraphics[width=0.45\\textwidth]{{{image}}}}}\n".format(caption=subcaption1, image=image1))
    outfp.write("\\subfloat[{caption}]{{\includegraphics[width=0.45\\textwidth]{{{image}}}}}\n".format(caption=subcaption2, image=image2))
    outfp.write("\\caption{{{caption}}}\n".format(caption=maincaption))
    outfp.write("\\end{figure}\n")
    outfp.write("\\end{quote}\n")

def Report_Images_Aboveeachother(image1, subcaption1, image2, subcaption2, maincaption, outfp):
    'Print 2 images side-by-side to the latex document'
    outfp.write("\\begin{figure}[!htpb]\n")
    outfp.write("\\begin{minipage}{\\textwidth}\n")
    outfp.write("\includegraphics[width=\\textwidth]{{{image}}}\n".format(image=image1))
    outfp.write("\\caption{{{caption}}}\n".format(caption=subcaption1))
    outfp.write("\\end{minipage}\n")
    outfp.write("\\begin{minipage}{\\textwidth}\n")
    outfp.write("\includegraphics[width=\\textwidth]{{{image}}}\n".format(image=image2))
    outfp.write("\\caption{{{caption}}}\n".format(caption=subcaption2))
    outfp.write("\\end{minipage}\n")
    outfp.write("\\end{figure}\n")

def Report_Table_Channels(outfp):
    'Print the table of channel statistics.'

    # Channel stats
    num_channels = 512
    num_active_channels = len(set(_E['channel_number']))
    pct_active_channels = '{0:,.1f}'.format(num_active_channels / float(num_channels) * 100.0)
    tmp1 = np.bincount(_E['channel_number'])
    tmp2 = np.nonzero(tmp1)[0]
    tmp3 = zip(tmp2,tmp1[tmp2])
    num_active_channels_that_produced_only_one_read = len([1 for x in tmp3 if x[1] == 1])
    channel_speed = '{0:,.1f}'.format(np.mean(_E['events_per_second']))

    # Pore stats
    num_pores = 2048
    elts = [((_E['channel_number'][i] * 10) +_E['start_mux'][i]) for i in range(0, len(_E['channel_number']))]
    num_active_pores = len(set(elts))
    pct_active_pores = '{0:,.1f}'.format(num_active_pores / float(num_pores) * 100.0)
    tmp1 = np.bincount(elts)
    tmp2 = np.nonzero(tmp1)[0]
    tmp3 = zip(tmp2,tmp1[tmp2])
    num_active_pores_that_produced_only_one_read = len([1 for x in tmp3 if x[1] == 1])
    pore_speed = '-'

    # Device stats
    num_events = sum(_E['event_count'])
    run_start = min(set(_E['exp_start_time'])) # min(_E['read_start_seconds'])
    run_finish =  run_start + max(_E['read_end_hours_since_expt_start'])*60*60 # max(_E['read_end_seconds'])
    total_seconds = max(_E['read_end_hours_since_expt_start'])*60*60 #run_finish - run_start
    device_speed = '{0:,.3f}'.format((num_events / float(total_seconds) * 60 * 60 / 1000000.0))

    outfp.write("\\begin{quote}\\\\\n")
    outfp.write("\\begin{table}[!htpb]\\\\\n")
    outfp.write("\\small\\\\\n")
    outfp.write("\\tabcolsep=0.2cm\\\\\n")
    outfp.write("\\centering\\\\\n")
    outfp.write("\\begin{tabular}{ R{1.9cm} R{1.6cm} R{1.6cm} R{1.4cm} R{1.5cm} R{1.6cm} }\\\\\n")
    outfp.write("\\hline\n")
    #outfp.write("channels \% active&channel num active&channel num total&1 read only&channel speed e/s&device speed Me/h\\\\\n")
    #outfp.write("\\hline\n")
    #outfp.write("{0}\\\\\n".format('&'.join([str(x) for x in [pct_active_channels, num_active_channels, num_channels, num_active_channels_that_produced_only_one_read, channel_speed, device_speed]])))

    R1 = ['channel', pct_active_channels, num_active_channels, num_channels,
          num_active_channels_that_produced_only_one_read, '{0} e/s'.format(channel_speed)]
    R2 = ['pores', pct_active_pores, num_active_pores, num_pores,
          num_active_pores_that_produced_only_one_read, pore_speed]
    R3 = ['device', '-', '-', '-', '-', '{0} Me/h'.format(device_speed)]
    outfp.write("&\% active&num active&num total&1 read only&speed\\\\\n")
    outfp.write("\\hline\n")
    outfp.write("{0}\\\\\n".format('&'.join([str(x) for x in R1])))
    outfp.write("{0}\\\\\n".format('&'.join([str(x) for x in R2])))
    outfp.write("{0}\\\\\n".format('&'.join([str(x) for x in R3])))
    outfp.write("\\hline\n")
    outfp.write("\\end{tabular}\\\\\n")
    outfp.write("\\end{table}\\\\\n")
    outfp.write("\\end{quote}\\\\\n")

    # Append more information to the _tabledata list.
    global _tabledata
    _tabledata.append( ['channelpctactive', pct_active_channels] )
    _tabledata.append( ['channelnumactive', num_active_channels] )
    _tabledata.append( ['channelnumtotal', num_channels] )
    _tabledata.append( ['channelsyielding1readonly', num_active_channels_that_produced_only_one_read] )
    _tabledata.append( ['channelspeedeventspersecond', channel_speed] )
    _tabledata.append( ['porepctactive', pct_active_pores] )
    _tabledata.append( ['porenumactive', num_active_pores] )
    _tabledata.append( ['porenumtotal', num_pores] )
    _tabledata.append( ['poresyielding1readonly', num_active_pores_that_produced_only_one_read] )
    _tabledata.append( ['devicespeedMeventsperhour', device_speed] )

def Report_Table_Yield(outfp):
    'Print the table of yield statistics.'

    min_events = [x for x in set(_C['min_events'])][0]
    max_events = [x for x in set(_C['max_events'])][0]

    Ecallmask = (_E['event_count'] >= min_events) & (_E['event_count'] <= max_events)
    twodpassmask = _C['twod_readclass'] == 'pass'
    _tabledata.append( ['', XXXX] )

    _tabledata.append( ['devicespeedMeventsperhour', device_speed] )

def Report_Table_Yield(outfp):
    'Print the table of yield statistics.'

    min_events = [x for x in set(_C['min_events'])][0]
    max_events = [x for x in set(_C['max_events'])][0]

    Ecallmask = (_E['event_count'] >= min_events) & (_E['event_count'] <= max_events)
    twodpassmask = _C['twod_readclass'] == 'pass'
    _tabledata.append( ['', XXXX] )

    _tabledata.append( ['devicespeedMeventsperhour', device_speed] )

def Report_Table_Yield(outfp):
    'Print the table of yield statistics.'

    min_events = [x for x in set(_C['min_events'])][0]
    max_events = [x for x in set(_C['max_events'])][0]

    Ecallmask = (_E['event_count'] >= min_events) & (_E['event_count'] <= max_events)
    twodpassmask = _C['twod_readclass'] == 'pass'
    _tabledata.append( ['', XXXX] )

    _tabledata.append( ['devicespeedMeventsperhour', device_speed] )

def Report_Table_Yield(outfp):
    'Print the table of yield statistics.'

    min_events = [x for x in set(_C['min_events'])][0]
    max_events = [x for x in set(_C['max_events'])][0]

    Ecallmask = (_E['event_count'] >= min_events) & (_E['event_count'] <= max_events)
    twodpassmask = _C['twod_readclass'] == 'pass'

    _tabledata.append( ['devicespeedMeventsperhour', device_speed] )

def Report_Table_Yield(outfp):
    'Print the table of yield statistics.'

    min_events = [x for x in set(_C['min_events'])][0]
    max_events = [x for x in set(_C['max_events'])][0]

    Ecallmask = (_E['event_count'] >= min_events) & (_E['event_count'] <= max_events)
    twodpassmask = _C['twod_readclass'] == 'pass'
    twodfailmask = _C['twod_readclass'] == 'fail'

    count_Ecall = len(_E[Ecallmask])
    count_Eall = len(_E)
    count_2Dpass = sum((_C['twod_fastq_seqlen'] > 0) & (_C['twod_readclass'] == 'pass'))
    count_2Dall = sum(_C['twod_fastq_seqlen'] > 0)
    count_2Dfail = count_2Dall - count_2Dpass
    count_comp = sum(_C['complement_fastq_seqlen'] > 0)
    count_temp = sum(_C['template_fastq_seqlen'] > 0)

    D = {}
    D['read count'] = {
        'E-call' : '{0:,d}'.format(count_Ecall),
        'E-all' : '{0:,d}'.format(count_Eall),
        '2D-pass' : '{0:,d}'.format(count_2Dpass),
        '2D-fail' : '{0:,d}'.format(count_2Dfail),
        '2D-all' : '{0:,d}'.format(count_2Dall),
        'comp' : '{0:,d}'.format(count_comp),
        'temp' : '{0:,d}'.format(count_temp)
    }
    D['yield (M)'] = {
        'E-call' : '{0:,.1f}'.format(sum(_E['event_count'][Ecallmask]) / 1000000.0),
        'E-all' : '{0:,.1f}'.format(sum(_E['event_count']) / 1000000.0),
        '2D-pass' : '{0:,.1f}'.format(sum(_C['twod_fastq_seqlen'][twodpassmask]) / 1000000.0),
        '2D-fail' : '{0:,.1f}'.format((sum(_C['twod_fastq_seqlen']) - sum(_C['twod_fastq_seqlen'][twodpassmask])) / 1000000.0),
        '2D-all' : '{0:,.1f}'.format(sum(_C['twod_fastq_seqlen']) / 1000000.0),
        'comp' : '{0:,.1f}'.format(sum(_C['complement_fastq_seqlen']) / 1000000.0),
        'temp' : '{0:,.1f}'.format(sum(_C['template_fastq_seqlen']) / 1000000.0)
    }
    D['% E-call'] = {
        'E-call' : '100',
        'E-all' : '.',
        '2D-pass' : '{0:,.1f}'.format(float(count_2Dpass) / count_Ecall * 100.0),
        '2D-fail' : '{0:,.1f}'.format(float(count_2Dfail) / count_Ecall * 100.0),
        '2D-all' : '{0:,.1f}'.format(float(count_2Dall) / count_Ecall * 100.0),
        'comp' : '{0:,.1f}'.format(float(count_comp) / count_Ecall * 100.0),
        'temp' : '{0:,.1f}'.format(float(count_temp) / count_Ecall * 100.0)
    }
    D['% temp'] = {
        'E-call' : '.',
        'E-all' : '.',
        '2D-pass' : '{0:,.1f}'.format(float(count_2Dpass) / count_temp * 100.0),
        '2D-fail' : '{0:,.1f}'.format(float(count_2Dfail) / count_temp * 100.0),
        '2D-all' : '{0:,.1f}'.format(float(count_2Dall) / count_temp * 100.0),
        'comp' : '{0:,.1f}'.format(float(count_comp) / count_temp * 100.0),
        'temp' : '100'
    }
    D['len mean'] = {
        'E-call' : '{0:,.0f}'.format(np.mean(_E['event_count'][(_E['event_count'] >= min_events) & (_E['event_count'] <= max_events) & (_E['event_count'] > 0)])),
        'E-all' : '{0:,.0f}'.format(np.mean(_E['event_count'][_C['event_count'] > 0])),
        '2D-pass' : '{0:,.0f}'.format(np.mean(_C['twod_fastq_seqlen'][(_C['twod_fastq_seqlen'] > 0) & (_C['twod_readclass'] == 'pass')]) \
                    if len(_C['twod_fastq_seqlen'][(_C['twod_fastq_seqlen'] > 0) & (_C['twod_readclass'] == 'pass')]) else 0),
        '2D-fail' : '{0:,.0f}'.format(np.mean(_C['twod_fastq_seqlen'][(_C['twod_fastq_seqlen'] > 0) & (_C['twod_readclass'] == 'fail')]) \
                    if len(_C['twod_fastq_seqlen'][(_C['twod_fastq_seqlen'] > 0) & (_C['twod_readclass'] == 'fail')]) else 0),
        '2D-all' : '{0:,.0f}'.format(np.mean(_C['twod_fastq_seqlen'][_C['twod_fastq_seqlen'] > 0]) \
                    if len(_C['twod_fastq_seqlen'][_C['twod_fastq_seqlen'] > 0]) else 0),
        'comp' : '{0:,.0f}'.format(np.mean(_C['complement_fastq_seqlen'][_C['complement_fastq_seqlen'] > 0]) \
                    if len(_C['complement_fastq_seqlen'][_C['complement_fastq_seqlen'] > 0]) else 0),
        'temp' : '{0:,.0f}'.format(np.mean(_C['template_fastq_seqlen'][_C['template_fastq_seqlen'] > 0]))
    }
    D['len mdn'] = {
        'E-call' : '{0:,.0f}'.format(np.median(_E['event_count'][(_E['event_count'] >= min_events) & (_E['event_count'] <= max_events) & (_E['event_count'] > 0)])),
        'E-all' : '{0:,.0f}'.format(np.median(_E['event_count'][_C['event_count'] > 0])),
        '2D-pass' : '{0:,.0f}'.format(np.median(_C['twod_fastq_seqlen'][(_C['twod_fastq_seqlen'] > 0) & (_C['twod_readclass'] == 'pass')]) \
                     if len(_C['twod_fastq_seqlen'][(_C['twod_fastq_seqlen'] > 0) & (_C['twod_readclass'] == 'pass')]) else 0),
        '2D-fail' : '{0:,.0f}'.format(np.median(_C['twod_fastq_seqlen'][(_C['twod_fastq_seqlen'] > 0) & (_C['twod_readclass'] == 'fail')]) \
                     if len(_C['twod_fastq_seqlen'][(_C['twod_fastq_seqlen'] > 0) & (_C['twod_readclass'] == 'fail')]) else 0),
        '2D-all' : '{0:,.0f}'.format(np.median(_C['twod_fastq_seqlen'][_C['twod_fastq_seqlen'] > 0]) \
                     if len(_C['twod_fastq_seqlen'][_C['twod_fastq_seqlen'] > 0]) else 0),
        'comp' : '{0:,.0f}'.format(np.median(_C['complement_fastq_seqlen'][_C['complement_fastq_seqlen'] > 0]) \
                     if len(_C['complement_fastq_seqlen'][_C['complement_fastq_seqlen'] > 0]) else 0),
        'temp' : '{0:,.0f}'.format(np.median(_C['template_fastq_seqlen'][_C['template_fastq_seqlen'] > 0]))
    }
    D['len min'] = {
        'E-call' : '{0:,.0f}'.format(min(_E['event_count'][(_E['event_count'] >= min_events) & (_E['event_count'] <= max_events) & (_E['event_count'] > 0)])),
        'E-all' : '{0:,.0f}'.format(min(_E['event_count'][_C['event_count'] > 0])),
        '2D-pass' : '{0:,.0f}'.format(min(_C['twod_fastq_seqlen'][(_C['twod_fastq_seqlen'] > 0) & (_C['twod_readclass'] == 'pass')]) \
                  if len(_C['twod_fastq_seqlen'][(_C['twod_fastq_seqlen'] > 0) & (_C['twod_readclass'] == 'pass')]) else 0),
        '2D-fail' : '{0:,.0f}'.format(min(_C['twod_fastq_seqlen'][(_C['twod_fastq_seqlen'] > 0) & (_C['twod_readclass'] == 'fail')]) \
                    if len(_C['twod_fastq_seqlen'][(_C['twod_fastq_seqlen'] > 0) & (_C['twod_readclass'] == 'fail')]) else 0),
        '2D-all' : '{0:,.0f}'.format(min(_C['twod_fastq_seqlen'][_C['twod_fastq_seqlen'] > 0]) if len(_C['twod_fastq_seqlen'][_C['twod_fastq_seqlen'] > 0]) else 0),
        'comp' : '{0:,.0f}'.format(min(_C['complement_fastq_seqlen'][_C['complement_fastq_seqlen'] > 0]) \
                  if len(_C['complement_fastq_seqlen'][_C['complement_fastq_seqlen'] > 0]) else 0),
        'temp' : '{0:,.0f}'.format(min(_C['template_fastq_seqlen'][_C['template_fastq_seqlen'] > 0]))
    }
    D['len max (K)' ] = {
        'E-call' : '{0:,.1f}'.format(max(_E['event_count'][(_E['event_count'] >= min_events) & (_E['event_count'] <= max_events) & (_E['event_count'] > 0)]) / 1000.0),
        'E-all' : '{0:,.1f}'.format(max(_E['event_count'][_C['event_count'] > 0]) / 1000.0),
        '2D-pass' : '{0:,.1f}'.format(max(_C['twod_fastq_seqlen'][(_C['twod_fastq_seqlen'] > 0) & (_C['twod_readclass'] == 'pass')]) / 1000.0 \
                    if len(_C['twod_fastq_seqlen'][(_C['twod_fastq_seqlen'] > 0) & (_C['twod_readclass'] == 'pass')]) else 0),
        '2D-fail' : '{0:,.1f}'.format(max(_C['twod_fastq_seqlen'][(_C['twod_fastq_seqlen'] > 0) & (_C['twod_readclass'] == 'fail')]) / 1000.0 \
                    if len(_C['twod_fastq_seqlen'][(_C['twod_fastq_seqlen'] > 0) & (_C['twod_readclass'] == 'fail')]) else 0),
        '2D-all' : '{0:,.1f}'.format(max(_C['twod_fastq_seqlen'][_C['twod_fastq_seqlen'] > 0]) / 1000.0 \
                    if len(_C['twod_fastq_seqlen'][_C['twod_fastq_seqlen'] > 0]) else 0),
        'comp' : '{0:,.1f}'.format(max(_C['complement_fastq_seqlen'][_C['complement_fastq_seqlen'] > 0]) / 1000.0 \
                    if len(_C['complement_fastq_seqlen'][_C['complement_fastq_seqlen'] > 0]) else 0),
        'temp' : '{0:,.1f}'.format(max(_C['template_fastq_seqlen'][_C['template_fastq_seqlen'] > 0]) / 1000.0)
    }
    D['bq mean'] = {
        'E-call' : '.',
        'E-all' : '.',
        '2D-pass' : '{0:,.1f}'.format(np.mean(_C['twod_fastq_bqmean'][(_C['twod_fastq_bqmean'] > 0) & (_C['twod_readclass'] == 'pass')]) \
                    if len(_C['twod_fastq_bqmean'][(_C['twod_fastq_bqmean'] > 0) & (_C['twod_readclass'] == 'pass')]) else 0),
        '2D-fail' : '{0:,.1f}'.format(np.mean(_C['twod_fastq_bqmean'][(_C['twod_fastq_bqmean'] > 0) & (_C['twod_readclass'] == 'fail')]) \
                    if len(_C['twod_fastq_bqmean'][(_C['twod_fastq_bqmean'] > 0) & (_C['twod_readclass'] == 'fail')]) else 0),
        '2D-all' : '{0:,.1f}'.format(np.mean(_C['twod_fastq_bqmean'][_C['twod_fastq_bqmean'] > 0]) \
                    if len(_C['twod_fastq_bqmean'][_C['twod_fastq_bqmean'] > 0]) else 0),
        'comp' : '{0:,.1f}'.format(np.mean(_C['complement_fastq_bqmean'][_C['complement_fastq_bqmean'] > 0]) \
                    if len(_C['complement_fastq_bqmean'][_C['complement_fastq_bqmean'] > 0]) else 0),
        'temp' : '{0:,.1f}'.format(np.mean(_C['template_fastq_bqmean'][_C['template_fastq_bqmean'] > 0])) 
    }
    D['bq mdn'] = {
        'E-call' : '.',
        'E-all' : '.',
        '2D-pass' : '{0:,.1f}'.format(np.median(_C['twod_fastq_bqmean'][(_C['twod_fastq_bqmean'] > 0) & (_C['twod_readclass'] == 'pass')]) \
                    if len(_C['twod_fastq_bqmean'][(_C['twod_fastq_bqmean'] > 0) & (_C['twod_readclass'] == 'pass')]) else 0),
        '2D-fail' : '{0:,.1f}'.format(np.median(_C['twod_fastq_bqmean'][(_C['twod_fastq_bqmean'] > 0) & (_C['twod_readclass'] ==  'fail')]) \
                    if len(_C['twod_fastq_bqmean'][(_C['twod_fastq_bqmean'] > 0) & (_C['twod_readclass'] == 'fail')]) else 0),
        '2D-all' : '{0:,.1f}'.format(np.median(_C['twod_fastq_bqmean'][_C['twod_fastq_bqmean'] > 0]) \
                    if len(_C['twod_fastq_bqmean'][_C['twod_fastq_bqmean'] > 0]) else 0),
        'comp' : '{0:,.1f}'.format(np.median(_C['complement_fastq_bqmean'][_C['complement_fastq_bqmean'] > 0]) \
                    if len(_C['complement_fastq_bqmean'][_C['complement_fastq_bqmean'] > 0]) else 0),
        'temp' : '{0:,.1f}'.format(np.median(_C['template_fastq_bqmean'][_C['template_fastq_bqmean'] > 0])) 
    }
    R = [ 'yield (M)', 'read count', '% E-call', '% temp', 'len mean', 'len mdn', 'len min', 'len max (K)', 'bq mean', 'bq mdn' ]
    H = [ 'E-call', 'E-all', '2D-pass', '2D-fail', '2D-all', 'comp', 'temp' ]

    outfp.write("\\begin{quote}\\\\\n")
    outfp.write("\\begin{table}[!htpb]\\\\\n")
    outfp.write("\\small\\\\\n")
    outfp.write("\\tabcolsep=0.2cm\\\\\n")
    outfp.write("\\centering\\\\\n")
    outfp.write("\\begin{tabular}{ L{1.7cm} R{1cm} R{1.2cm} R{1.2cm} R{1cm} R{1cm} R{1cm} R{1.5cm} }\\\\\n")
    outfp.write("\\hline\n")
    outfp.write("{0}\\\\\n".format('&'.join([''] + H)))
    outfp.write("\\hline\n")
    for rowname in R:
        outfp.write("{0}\\\\\n".format('&'.join([escaped(rowname)] + [D[rowname][x] for x in H])))
    outfp.write("\\hline\n")
    outfp.write("\\end{tabular}\\\\\n")
    outfp.write("\\end{table}\\\\\n")
    outfp.write("\\end{quote}\\\\\n")

    # Append the yield information, one cell per line.
    global _tabledata
    for readtype in H:
        for stattype in R:
            _tabledata.append( ['{0} {1}'.format(readtype, stattype), D[stattype][readtype]] )

def Report_DataSummary():
    _outfp['report'].write("\\section{Data summary}\n\n")

    _outfp['report'].write("\\subsection{Experiment summary}\n\n")
    Report_Table(sectionD['summary'], "Summary", True, _outfp['report'])
    _outfp['report'].write("\\subsection{Sample preparation, run configuration and bioinformatics QC parameters}\n\n")
    Report_Table(sectionD['sample'], "Sample preparation", True, _outfp['report'])
    _outfp['report'].write("\\subsection{Flowcell performance and data yield}\n\n")
    Report_Table_Channels(_outfp['report'])
    Report_Table_Yield(_outfp['report'])
    #_outfp['report'].write("\\subsection{Quantity of target genome(s) acquired}\n\n")

def Report_ReadPlots():

    _outfp['report'].write("\\raggedright\n")
    _outfp['report'].write("\\section{Read summary}\n\n")

    _outfp['report'].write("\\raggedright\n")
    _outfp['report'].write("\\subsection{Sequencing statistics}\n\n")
    pngdetails = [
        [ 'runtime_readcompletioncount', 'Read yield', '' ],
        [ 'runtime_readcompletionrate', 'Read completion rate', '' ],
        [ 'runtime_eventcompletioncount', 'Event yield', '' ],
        [ 'runtime_eventcompletionrate', 'Event completion rate', '' ],
        [ 'runtime_activechannels', 'Channel activity', '' ],
       #[ 'runtime_channelspeedratio', 'Channel speed ratio', 'A boxplot of the ratio of the channel speed for the first and second half of all completed reads, where using "half" is an approximation for where the template and complement parts of the read might be, since we do not have this information before the read is basecalled.' ],
        [ 'runtime_channelspeedovertime', 'Channel speed over time', '' ],
        [ 'runtime_2Dmeanbqovertime', 'Mean BQ of 2D reads', '' ],
       #[ 'runtime_tempmeanbqovertime', 'Mean BQ of template reads', 'Ideally, this would be constant across a run.' ],
       #[ 'runtime_compmeanbqovertime', 'Mean BQ of complement reads', 'Ideally, this would be constant across a run.' ],
        [ 'runtime_chiptemp', 'Flow-cell temperature', '' ],
       #[ 'runtime_pctNin2Dreadsovertime', 'Percent of uncalled bases in 2D reads', 'Ideally, this would be zero.' ]
    ]
    for image in pngdetails:
        pngpath = os.path.join('pngpdf', '{0}_{1}.png'.format(_M['outfileprefix'], image[0]))
        Report_Image(pngpath, image[1], image[2], _outfp['report'])

    _outfp['report'].write("\\raggedright\n")
    _outfp['report'].write("\\subsection{Read length statistics}\n\n")
    pngdetails = [
        [ 'readlen_histlower', 'Length of event-reads around the median', 'Distribution of most of the event-read lengths.' ]
    ]
    for image in pngdetails:
        pngpath = os.path.join('pngpdf', '{0}_{1}.png'.format(_M['outfileprefix'], image[0]))
        Report_Image(pngpath, image[1], image[2], _outfp['report'])
    pngdetails = [
        [ 'readlen_histall', 'Length of all event-reads', 'Distribution of all event-read lengths generated during the sequencing run.' ],
        [ 'readlen_histupper', 'Length of upper trailing event-read lengths', 'Distribution of the higher event-read lengths.' ]
    ]
    image1 = os.path.join('pngpdf', '{0}_{1}.png'.format(_M['outfileprefix'], pngdetails[0][0]))
    image2 = os.path.join('pngpdf', '{0}_{1}.png'.format(_M['outfileprefix'], pngdetails[1][0]))
    subcaption1 = 'All event-reads' # pngdetails[0][1]
    subcaption2 = 'Trailing event-reads only' # pngdetails[1][1]
    maincaption = 'Event-read length distributions'
    Report_Images_Sidebyside(image1, subcaption1, image2, subcaption2, maincaption, _outfp['report'])
    pngdetails = [
        [ 'readlen_yield', 'Event yield for all event-reads', 'MinKNOW-like bar-plot showing the number of events in reads of different lengths.' ],
        [ 'readlen_yieldlower', 'Event yield for most of the event-reads', 'MinKNOW-like bar-plot showing the number of events in reads of different lengths, excluding the right tail.' ],
       #[ 'readlen_channelspeed', 'Channel speed along the reads', '10 boxplots showing the channel speed for each equal-sized segment of read. Although a rough estimate, this might take into account a difference in channel speed before and after the hairpin , which is expected to be near the middle of the read.' ],
       #[ 'readlen_durationvseventcount', 'Channel speed for reads of different lengths', 'Scatter plot that shows the relationship between how long the read took to generate and how much data was generated, that is, if shorter/longer reads go through the pore faster/slower.' ]
    ]
    for image in pngdetails:
        pngpath = os.path.join('pngpdf', '{0}_{1}.png'.format(_M['outfileprefix'], image[0]))
        Report_Image(pngpath, image[1], image[2], _outfp['report'])

    _outfp['report'].write("\\raggedright\n")
    _outfp['report'].write("\\subsection{Channel performance statistics}\n\n")
    pngdetails = [
        [ 'channel_readcountcallable', 'Channel yield of event-reads with a callable length', '' ],
        [ 'channel_readcounttooshort', 'Channel yield of event-reads too short to call', '' ],
        [ 'channel_readcounttoolong', 'Channel yield of event-reads too long to call', '' ],
       #[ 'channel_readlenboxplot', 'Read lengths produced by each channel', 'Boxplot of the read lengths produced by each active channel.' ],
        [ 'channel_readlenmean', 'Mean read lengths produced by each channel', '' ],
        [ 'channel_initfilldelay', 'Initial fill delay of each channel', '' ],
        [ 'channel_refilldelay', 'Re-fill delay of each channel', '' ],
        [ 'channel_activitytime', 'Channel activity over time', 'Colours denote the length of the resulting event-read.' ],
    ]
    for image in pngdetails:
        pngpath = os.path.join('pngpdf', '{0}_{1}.png'.format(_M['outfileprefix'], image[0]))
        Report_Image(pngpath, image[1], image[2], _outfp['report'])

    _outfp['report'].write("\\raggedright\n")
    _outfp['report'].write("\\subsection{Base-calling statistics}\n\n")
    pngdetails = [
        [ 'basecall_baseyield', 'Read data lost at each step', '' ],
        [ 'basecall_readyield', 'Read count lost at each step', '' ]
    ]
    image1 = os.path.join('pngpdf', '{0}_{1}.png'.format(_M['outfileprefix'], pngdetails[0][0]))
    image2 = os.path.join('pngpdf', '{0}_{1}.png'.format(_M['outfileprefix'], pngdetails[1][0]))
    subcaption1 = '' # pngdetails[0][1]
    subcaption2 = '' # pngdetails[1][1]
    maincaption = 'Read yield and count lost at each step'
    Report_Images_Sidebyside(image1, subcaption1, image2, subcaption2, maincaption, _outfp['report'])
    pngdetails = [
       #[ 'basecall_distribeventcountgivingnoTC2D', 'Event lengths that gave no 2D, complement or template reads', '' ],
        [ 'basecall_readlendistall', 'Distribution of called read lengths', '' ],
        [ 'basecall_readlendist2D', 'Distribution of 2D read lengths', '' ],
        [ 'basecall_readlendisttemp', 'Distribution of template read lengths', '' ],
        [ 'basecall_readlendistcomp', 'Distribution of complement read lengths', '' ]
    ]
    for image in pngdetails:
        pngpath = os.path.join('pngpdf', '{0}_{1}.png'.format(_M['outfileprefix'], image[0]))
        Report_Image(pngpath, image[1], image[2], _outfp['report'])
    pngdetails = [
        [ 'basecall_channelspeedtemp2comp', 'Channel speed of template and complement', '' ],
        [ 'basecall_eventto2Dreadlen', 'Relationship between complement and 2D read length', '' ]
    ]
    image1 = os.path.join('pngpdf', '{0}_{1}.png'.format(_M['outfileprefix'], pngdetails[0][0]))
    image2 = os.path.join('pngpdf', '{0}_{1}.png'.format(_M['outfileprefix'], pngdetails[1][0]))
    subcaption1 = '' # pngdetails[0][1]
    subcaption2 = '' # pngdetails[1][1]
    maincaption = 'Relationship between event and basecall-read lengths'
    Report_Images_Sidebyside(image1, subcaption1, image2, subcaption2, maincaption, _outfp['report'])
       #[ 'basecall_eventtotempreadlen', 'Relationship between event count and template read length', '' ],
       #[ 'basecall_eventtocompreadlen', 'Relationship between event count and complement read length', '' ],
       #[ 'basecall_eventtotempcompreadlen', 'Relationship between event count and template+complement read length', '' ],
       #[ 'basecall_tempto2Dreadlen', 'Relationship between template and 2D read lengths', '' ],
       #[ 'basecall_compto2Dreadlen', 'Relationship between event count and 2D read length', '' ],
       #[ 'basecall_tempcompto2Dreadlen', 'Relationship between read length of temp+comp and 2D read', '' ]


def Report_RunPerformance():
    ''
    plotL = [
        # How many molecules read and how long were they?
        'Sequencing success rate of each channel (interleaved bar-plot of number of successful and unsuccessful molecule read attempts per channel)',
        'Number of ssDNA molecules produced by each channel (bar-plot of number of event files for each channel, report number and percent of active and dead channels)',
        'Length of ssDNA molecules sequenced successfully (histogram of event lengths across run, and report median)',
        'Length of ssDNA molecules sequenced successfully, without tail (histogram or event lengths across run, up to median*2)',
        'Time that each ssDNA molecule was read (line plots, one per channel, xaxis is time, yaxis is channel number, horizontal lines to show when template and complement were being read - use different colours for the T and C read components)',

        'Number of template, complement and 2D reads from channel (interleaved bar-plot or 3 line-plots)',
        'Relationship between number of events and number of bases yielded (scatterplot of events on xaxis and temp+comp bases on yaxis, and fit a linear regression line through it, report slope and intercept formula)',

        'Distribution of template, complement and 2D read lengths (3 histograms, report median of each)',
        'Data yield during the run (bar-plot of total number of events produced at regular intervals during run)',
        '',
        'Mean base quality along read (bar-plot for each base, try to have line plot of number of reads on which boxplot is based on same plot and its own axis on RHS) - or insert Per Base Sequence Quality from FastQC',
        'Mean sequence quality (quality score distribution ove all sequences - Per Sequence Quality Scores from FastQC?)',
        'Per base sequence content (shows proportion of ACGT against position in read - FastQC)',
        'per base GC content (shows GC bias for position in read - FastQC)',
        'Per sequence GC content - GC distrib over all sequences - FastQC)',
        'Per base N content - Position of Ns in a read - FastQC',
        '',
        ''
    ]
    _outfp['report'].write("\\section{Run performance}\n\n")
    Report_Table(_D['runperformance'], "Sequencing run statistics.", True, _outfp['report'])
    Report_BulletList(plotL, _outfp['report'])
    _outfp['report'].write('\\nopagebreak[2]\n')

def Report_DataYield():
    ''
    _outfp['report'].write("\\section{Data yield}\n\n")
    Report_Table(_D['datayield'], "Data yield statistics.", True, _outfp['report'])
    _outfp['report'].write('\\nopagebreak[2]\n')

def Report_DataQuality():
    ''
    _outfp['report'].write("\\section{Data quality}\n\n")
    Report_Table(_D['dataquality'], "Data quality statistics.", True, _outfp['report'])
    _outfp['report'].write('\\nopagebreak[2]\n')

def Report_TargetStats():
    ''
    _outfp['report'].write("\\section{Target statistics based on mapping to a reference}\n\n")
    Report_Table(_D['targetstats'], "Target genome statistics.", True, _outfp['report'])
    _outfp['report'].write('\\nopagebreak[2]\n')

def Report_Close():
    'Report closing statements.'
    global _outfp
    _outfp['report'].write("\end{document}\n")
    _outfp['report'].close()

def Report_ConvertToPDF():
    'Convert the .tex to .pdf.'
    cmd = "tex2pdf {texfile} &>/dev/null".format(texfile=_outpath['report'])
    os.system(cmd)

def Create_Report():
    'Create a PDF report using LaTeX.'

    if os.path.exists(_outpath['reportpdf']) and not  _args.overwrite and not _args.forcereportstep:
        sys.stdout.write('Info: PDF report already up to date, no need to regenerate\n')
        return
    sys.stdout.write('Info: PDF report will be regenerated\n')

    Report_Init()
    Report_DataSummary()
    Report_ReadPlots()
    Report_Close()
    Report_ConvertToPDF()

# ============================================================================ #
# Save all data required for new reports, one file per table or plot.          #
# ============================================================================ #

def Speed_Desc(n):
    if n < 20:
        desc = 'POOR <20'
    elif n < 30:
        desc = 'AVERAGE >=20 but <30'
    else:
        desc = 'GOOD >=30'
    return desc

def Save_TableData():
    'Save the data that will be printed in tabular format at the start of the report.'

    global sectionD

    if os.path.exists(_outpath['reportpdf']) and not  _args.overwrite and not _args.forcereportstep:
        sys.stdout.write('Info: Table data for .tex report already up to date\n')
        return
    sys.stdout.write('Info: Table data for .tex report being updated\n')

    sectionT = [
        ['summary',  'PoreQC sequencing report'],
        ['sample', 'Sample preparation and run configuration'],
        #['sequencing', 'Sequencing run summary'],
        ['basecalling', 'Base-calling summary']
        #['mapping', 'Quantity of target genome(s) acquired']
    ]

    #run_start = min(_E['read_start_seconds'])
    #run_finish = max(_E['read_end_seconds'])
    run_start = min(set(_E['exp_start_time'])) # min(_E['read_start_seconds'])
    run_finish =  run_start + max(_E['read_end_hours_since_expt_start'])*60*60 # max(_E['read_end_seconds'])
    run_duration = datetime.timedelta(seconds=run_finish) - datetime.timedelta(seconds=run_start)

    weeks, days = divmod(run_duration.days, 7)
    minutes, seconds = divmod(run_duration.seconds, 60)
    hours, minutes = divmod(minutes, 60)
    run_duration_text = "{d} day{pl} {h:02d}:{m:02d}:{s:02d}".format(d=days, pl='s' if days>1 else '', h=hours, m=minutes, s=seconds)

    report_type = 'real-time progress report' if _S == 'Running' else 'final sequencing run report'
    sectionD = {}
    sectionD['summary'] = []
    sectionD['summary'].append(['Experiment aims', _M['experimentaims']])
    sectionD['summary'].append(['', ''])
    sectionD['summary'].append(['Experiment description', _M['experimentdesc']])
    sectionD['summary'].append(['', ''])
    sectionD['summary'].append(['Sample name (short)', _M['samplenameshort']])
    sectionD['summary'].append(['Output file prefix', _M['outfileprefix']])
    sectionD['summary'].append(['', ''])
    sectionD['summary'].append(['Experiment start', time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(run_start))])
    sectionD['summary'].append(['Experiment finish', time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(run_finish))])
    sectionD['summary'].append(['Experiment duration', run_duration_text])
    sectionD['summary'].append(['', ''])

    workflow_name = ','.join([x.replace('.py', '') for x in set(_C['workflow_script'])])
    version = ','.join([x for x in set(_C['nameversion'])])
    workflow_fullname = ','.join([x for x in set(_C['workflow_fullname'])])
    workflow_shortname = ','.join([x for x in set(_C['workflow_shortname'])])
    workflow_version = ','.join([x for x in set(_C['workflow_version'])])
    model_type = ','.join([x for x in set(_C['model_type'])])
    metrichorworkflows_version = '{0} {1} {2}'.format(workflow_fullname, workflow_shortname, model_type).strip(',').strip(' ')
    chemistry_version_metadata = _M['chemistryversion'].strip()
    try:
        chemistry_version_files = ','.join([x for x in set(_C['model_type'])])
    except:
        chemistry_version_files = ''
    chemistry_version = chemistry_version_metadata
    if len(chemistry_version_files) and chemistry_version_metadata.lower() != chemistry_version_files.lower() and chemistry_version_files.lower() != 'unknown':
       chemistry_version += ' ({0} in files)'.format(chemistry_version_files)
    min_eventsS = ','.join(['{0:,}'.format(x) for x in set(_C['min_events'])])
    max_eventsS = ','.join(['{0:,}'.format(x) for x in set(_C['max_events'])])
    num_events = sum(_E['event_count'])

    sectionD['sample'] = []
    sectionD['sample'].append(['Sample name (full)', _M['samplenamefull']])
    sectionD['sample'].append(['Species', _M['species']])
    sectionD['sample'].append(['Flowcell version', chemistry_version])
    sectionD['sample'].append(['', ''])
    sectionD['sample'].append(['Flowcell ID', _E[0]['flow_cell_id'] if len(_E[0]['flow_cell_id']) else _M['flowcellid']])
    sectionD['sample'].append(['Device ID', _E[0]['device_id']])
    sectionD['sample'].append(['Program length (hrs)', _M['programlength_hrs']])
    sectionD['sample'].append(['MinKNOW version', _M['minknowversion']])
    sectionD['sample'].append(['Metrichor version', _metrichor_version])
    sectionD['sample'].append(['Base-calling workflow', metrichorworkflows_version])
    sectionD['sample'].append(['Base-callable event lengths', '{0} to {1}'.format(min_eventsS, max_eventsS)])
    sectionD['sample'].append(['', ''])
#    if not len(_R):
#        sectionD['sample'].append(['Reference for QC', 'None specified'])
#        #sectionD['sample'].append(['Reference length (bp)', 'Not applicable'])
#    else:
#        ref_size = sum([len(_R[id].seq) for id in _R.keys()])
#        theoretical_coverage = num_events / float(ref_size) * 0.5
#        sectionD['sample'].append(['Reference for QC', '{0} contig(s): {1}'.format(
#            len(_R.keys()), ';'.join([_R[id].description.replace('\t', ' ') for id in _R.keys()]))])
#        sectionD['sample'].append(['Reference length (bp)', '{0:,}'.format(sum([len(_R[id].seq) for id in _R.keys()]))])
#        sectionD['sample'].append(['Target reference GC%', gc_percent])
#        sectionD['sample'].append(['Theoretical target coverage', '{0:,.1f}-fold'.format(theoretical_coverage)])

    # Some bits are commented out because I haven't got the part based on basecalls working yet.
    num_channels = 512
    num_active_channels = len(set(_E['channel_number']))
    pct_active_channels = num_active_channels / float(num_channels) * 100.0
    tmp1 = np.bincount(_E['channel_number'])
    tmp2 = np.nonzero(tmp1)[0]
    tmp3 = zip(tmp2,tmp1[tmp2])
    num_active_channels_that_produced_only_one_read = len([1 for x in tmp3 if x[1] == 1])
    total_seconds = run_finish - run_start
    channel_speed = np.mean(_E['events_per_second'])
    channel_speed_desc = Speed_Desc(channel_speed)
    device_speed = (num_events / float(total_seconds) * 60 * 60 / 1000000.0)
    inreadsdir = os.path.join(os.path.expandvars(_args.inrundir), 'reads')
    reads_completed = len(_E)
    reads_inprogress = len([x for x in os.listdir(inreadsdir) if os.path.isfile(os.path.join(inreadsdir, x)) and x.endswith('.fast5.tmp')])
    if not len(_R):
        gc_percent = 'Not applicable'
        ref_size = 'Not applicable'
        #theoretical_coverage = 'Not applicable'
    else:
        gc_percent = '{0:.1f}'.format(GC(''.join([str(_R[id].seq) for id in _R.keys()])))
        ref_size = sum([len(_R[id].seq) for id in _R.keys()])
        #theoretical_coverage = num_events / float(ref_size) * 0.5
        #min_events = max(_C['min_events'])
        #max_events = min(_C['max_events'])
        #sectionD['sequencing'] = []
        #sectionD['sequencing'].append(['Target reference size', '{0:,}'.format(ref_size)])
        #sectionD['sequencing'].append(['Target reference GC%', gc_percent])
        #sectionD['sequencing'].append(['Theoretical target coverage', '{0:,.1f}-fold'.format(theoretical_coverage)])

    if not len(_R):
        sectionD['sample'].append(['Reference for QC', 'None specified'])
        #sectionD['sample'].append(['Reference length (bp)', 'Not applicable'])
    else:
        ref_size = sum([len(_R[id].seq) for id in _R.keys()])
        theoretical_coverage = num_events / float(ref_size) * 0.5
        sectionD['sample'].append(['Reference for QC', '{0} contig(s): {1}'.format(
            len(_R.keys()), ';'.join([_R[id].description.replace('\t', ' ') for id in _R.keys()]))])
        sectionD['sample'].append(['Reference length (bp)', '{0:,}'.format(sum([len(_R[id].seq) for id in _R.keys()]))])
        sectionD['sample'].append(['Target reference GC%', gc_percent])
        sectionD['sample'].append(['Theoretical target coverage', '{0:,.1f}-fold'.format(theoretical_coverage)])

    calls_completed = len(_C)
    num_template = sum(_C['template_fastq_seqlen'] > 0)
    pct_eventreads = float(num_template) / calls_completed * 100
    num_complement = sum(_C['complement_fastq_seqlen'] > 0)
    num_2D = sum(_C['twod_fastq_seqlen'] > 0)
    pct_complement = num_complement / float(num_template) * 100.0
    pct_2D = num_2D / float(num_template) * 100.0
    sectionD['basecalling'] = []
    sectionD['basecalling'].append(['Reads with base-call data', '{0:,} of {1:,} ({2:.1f}%)'.format(calls_completed, reads_completed, calls_completed / float(reads_completed) * 100.0, )])
    sectionD['basecalling'].append(['template', '{0:,} ({1:.1f}% of called event reads)'.format(num_template, pct_eventreads)])
    sectionD['basecalling'].append(['complement', '{0:,} ({1:.1f}% of template)'.format(num_complement, pct_complement)])
    sectionD['basecalling'].append(['2D', '{0:,} ({1:,.1f}% of template)'.format(num_2D, pct_2D)])
    sectionD['basecalling'].append(['Longest calls', ''])
    sectionD['basecalling'].append(['template', '{0:,}'.format(max(_C['template_fastq_seqlen']))])
    sectionD['basecalling'].append(['complement', '{0:,}'.format(max(_C['complement_fastq_seqlen']))])
    sectionD['basecalling'].append(['2D', '{0:,}'.format(max(_C['twod_fastq_seqlen']))])
    sectionD['basecalling'].append(['Mean,Median base quality', ''])
    L = _C['template_fastq_bqmean'][_C['template_fastq_bqmean'] != 0]
    sectionD['basecalling'].append(['template', '{0:.1f} , {1:.1f}'.format(
        np.mean(L) if len(L) else 0.0, np.median(L) if len(L) else 0.0)])
    L = _C['complement_fastq_bqmean'][_C['complement_fastq_bqmean'] != 0]
    sectionD['basecalling'].append(['complement', '{0:.1f} , {1:.1f}'.format(
        np.mean(L) if len(L) else 0.0, np.median(L) if len(L) else 0.0)])
    L = _C['twod_fastq_bqmean'][_C['twod_fastq_bqmean'] != 0]
    sectionD['basecalling'].append(['2D', '{0:.1f} , {1:.1f}'.format(
        np.mean(L) if len(L) else 0.0, np.median(L) if len(L) else 0.0)])
    sectionD['basecalling'].append(['Base-calling workflow', metrichorworkflows_version])
    sectionD['basecalling'].append(['Base-callable event lengths', '{0} to {1}'.format(min_eventsS, max_eventsS)])

    #sectionD['mapping'] = []
    #sectionD['mapping'].append(['Target genome length', '{0:,}'.format(ref_size)])
    #sectionD['mapping'].append(['Target genome GC%', gc_percent])
    #sectionD['mapping'].append(['Theoretical target coverage', '{0:,.1f}-fold'.format(theoretical_coverage)])
    #sectionD['mapping'].append(['Min read depth required', _M['mindepthrequired']])
    #sectionD['mapping'].append(['% genome > mindepth', ''])
    #sectionD['mapping'].append(['Mean template depth', ''])

    # Start accumulating the data that will be saved to two tabular files,
    # a .txt file with two columns that is human readable and a .csv file with
    # a header row and a data row that could be concatenated from the files of other
    # rows and loaded to a database or Excel for further processing later.

    global _tabledata
    _tabledata = []
    _tabledata.append( ['runfolder', os.path.basename(_args.inrundir)] )
    _tabledata.append( ['outprefix', _M['outfileprefix']] )
    _tabledata.append( ['runname', _M['samplenameshort']] )
    _tabledata.append( ['reportversion', _version] )
    _tabledata.append( ['reportrundt', time.strftime('%Y-%m-%d %H:%M:%S', _currentdatetime)] )
    _tabledata.append( ['exptname', _M['experimentname']] )
    _tabledata.append( ['exptstartdt', time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(run_start))] )
    _tabledata.append( ['exptfinishdt', time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(run_finish))] )
    _tabledata.append( ['exptdurationseconds', run_finish - run_start] )
    _tabledata.append( ['exptdurationtext', run_duration_text] )
    _tabledata.append( ['samplenamefull', _M['samplenamefull']] )
    _tabledata.append( ['samplespecies', _M['species']] )
    _tabledata.append( ['flowcellversion', chemistry_version] )
    _tabledata.append( ['flowcellid', _E[0]['flow_cell_id'] if len(_E[0]['flow_cell_id']) else _M['flowcellid']] )
    _tabledata.append( ['deviceid', _E[0]['device_id']] )
    _tabledata.append( ['programlenhrs', _M['programlength_hrs']] )
    _tabledata.append( ['minknowversion', _M['minknowversion']] )
    _tabledata.append( ['metrichorversion', _metrichor_version] )
    _tabledata.append( ['workflowtext', metrichorworkflows_version] )
    _tabledata.append( ['workflowfullname', _C[0]['workflow_fullname']] )
    _tabledata.append( ['workflowshortname', _C[0]['workflow_shortname']] )
    _tabledata.append( ['workflowversion', _C[0]['workflow_version']] )
    _tabledata.append( ['workflowmodeltype', model_type] )
    _tabledata.append( ['bcallminlen', min_eventsS] )
    _tabledata.append( ['bcallmaxlen', max_eventsS] )

def Save_ReadLengthData():
    '''
    Save the two tables of data needed to plot histogram of read lengths and
    total number of events in read length buckets (like displayed in MinKNOW).
    '''
    # Dump all the read lengths to a file.
    with open(_outpath['readlen'], 'w') as readlength_fp:
        readlength_fp.write('event_count\n')
        for val in _E[_E['event_count'] != 0]['event_count']:
            readlength_fp.write('{0}\n'.format(val))
    # Accumulate the read lengths into buckets of size --readlenbucketsz
    readlen_max = max(_E[_E['event_count'] != 0]['event_count'])
    readlen_bucketmax = int(math.ceil(readlen_max / float(_args.readlenbucketsz))) * _args.readlenbucketsz
    num_buckets = readlen_bucketmax / _args.readlenbucketsz
    eventcountbyreadlen = np.zeros((num_buckets,), dtype=_readlenaccum_dtype)
    for bucketidx in range(0, num_buckets):
        eventcountbyreadlen[bucketidx]['readlenmx'] = (bucketidx+1) * _args.readlenbucketsz
    for event_count in _E[_E['event_count'] != 0]['event_count']:
        bucketidx = int(math.floor(event_count / float(_args.readlenbucketsz)))
        event_count_in_mbases = float(event_count) / 1000000.0
        eventcountbyreadlen[bucketidx]['eventtotal_mbases'] = eventcountbyreadlen[bucketidx]['eventtotal_mbases'] + event_count_in_mbases
    with open(_outpath['readlenyield'], 'w') as readlenaccum_fp:
        readlenaccum_fp.write('readlength_max\ttotal_event_count\n')
        for bucketidx in range(0, len(eventcountbyreadlen)):
            readlenaccum_fp.write('{0}\t{1}\n'.format(eventcountbyreadlen[bucketidx]['readlenmx'], eventcountbyreadlen[bucketidx]['eventtotal_mbases']))

def Save_ReadStatData():
    '''
    Save a table listing the ch, file, channel_number, read_number, event_count, length_template, length_complement, length_2D
    for all the molecules that have been base-called.
    '''
    with open(_outpath['readstats'], 'w') as out_fp:
        out_fp.write('{0}\n'.format('\t'.join([x[0] for x in _readstats_dtype])))
        for rec in _B:
            if rec['read_length_template'] > 0:
                read_start_time = rec['exp_start_time'] + rec['start_time_read']
                read_end_time = read_start_time + rec['duration_read']
                template_start_time = rec['exp_start_time'] + rec['start_time_template']
                template_end_time = template_start_time + rec['duration_template']
                complement_start_time = rec['exp_start_time'] + rec['start_time_complement']
                complement_end_time = complement_start_time + rec['duration_complement']
                L = [
                    rec['ch'], rec['file'], rec['channel_number'], rec['read_number'], rec['event_count'],
                    rec['read_length_template'], rec['read_length_complement'], rec['read_length_2D'], rec['duration_read'], rec['duration_template'],
                    rec['duration_complement'], read_start_time, read_end_time, template_start_time, template_end_time,
                    complement_start_time, complement_end_time]
                out_fp.write('{0}\n'.format('\t'.join([str(x) for x in L])))


def Save_ChannelStatData():
    'Save this table, one row per channel, to a file.'

    channelstats = np.zeros((512+1,), _channelstats_dtype)
    for i in range(1, len(channelstats)):
         channelstats[i]['channel_number'] = i

    for rec in _B:
        channel_number = rec['channel_number']
        channelstats[channel_number]['num_reads'] += 1
        channelstats[channel_number]['num_template_reads'] += rec['has_template_read']
        channelstats[channel_number]['num_complement_reads'] += rec['has_complement_read']
        channelstats[channel_number]['num_2D_reads'] += rec['has_2D_read']
        channelstats[channel_number]['sum_event_count_mbases'] = channelstats[channel_number]['sum_event_count_mbases'] + rec['event_count']/1000000.0
        channelstats[channel_number]['sum_template_lengths_mbases'] = channelstats[channel_number]['sum_template_lengths_mbases'] + rec['read_length_template']/1000000.0
        channelstats[channel_number]['sum_complement_lengths_mbases'] = channelstats[channel_number]['sum_complement_lengths_mbases'] + rec['read_length_complement']/1000000.0
        channelstats[channel_number]['sum_2D_lengths_mbases'] = channelstats[channel_number]['sum_2D_lengths_mbases'] + rec['read_length_2D']/1000000.0

    with open(_outpath['channelstats'], 'w') as out_fp:
        out_fp.write('{0}\n'.format('\t'.join([x[0] for x in _channelstats_dtype])))
        for i in range(1, len(channelstats)):
            out_fp.write('{0}\n'.format('\t'.join([str(x) for x in channelstats[i]])))

def Save_RunStatusData():
    '''
    Collate the data required for doing time-series plots.
    For most of the plots, need a table which has one row per read molecule,
    the values of the metrics of interest, and the bucket number (category)
    to which is belongs. And then we can let R do the rest of the messy stuff.
    '''
    #timestats = np.zeros((N,), _timestats_dtype)

def Create_TabularDataFiles():
    if os.path.exists(_outpath['tabletxt']) and \
    os.path.exists(_outpath['tablecsv']) and \
    not _args.overwrite and not _args.forcereportstep:
        sys.stdout.write('Info: .txt and .csv tables already up to date, no need to regenerate\n')
        return
    sys.stdout.write('Info: .txt and .csv tables will be generated\n')

    with open(_outpath['tabletxt'], 'w') as out_fp:
        for elt in _tabledata:
            out_fp.write('{0}\t{1}\n'.format(elt[0], str(elt[1]).replace(',', '')))

    with open(_outpath['tablecsv'], 'w') as out_fp:
        out_fp.write('{0}\n'.format('\t'.join([x[0] for x in _tabledata])))
        out_fp.write('{0}\n'.format('\t'.join([str(x[1]).replace(',', '') for x in _tabledata])))

# ============================================================================ #
# Main                                                                         #
# ============================================================================ #

if __name__ == '__main__':

    Initialise()
    Read_Status()
    Read_Metadata()
    Read_Reference()
    Extract_NewEventData()
    Read_AllEventData()
    Extract_NewCalldata()
    Read_AllCallData()
    Save_TableData()
    Create_Plots()
    Create_Report()
    Create_TabularDataFiles()

# ============================================================================ #
