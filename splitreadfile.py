#!/usr/bin/env python

import numpy as np, os, sys
from Bio import SeqIO

if len(sys.argv) < 5:
    sys.stdout.write('Usage: splitreadfile.py in2Dallpath incallstatspath outdir outprefix\n'.format())
    sys.exit(1)
in2Dallpath = os.path.expandvars(sys.argv[1])
incallstatspath = os.path.expandvars(sys.argv[2])
outdir = os.path.expandvars(sys.argv[3])
outprefix = sys.argv[4]

linewidth = 100

if __name__ == '__main__':

    # Infer file type of in2Dallpath
    if in2Dallpath.endswith('.fasta'):
        infiletype = 'fasta'
    elif in2Dallpath.endswith('.fastq'):
        infiletype = 'fastq'
    else:
        sys.stderr.write('Erro: File of unknown suffix type ({0})\n'.format(in2Dallpath))
        sys.exit(1)

    # Read in the callstats.txt file (insplitinfopath) and work out which elements are high quality or not,
    # where high-quality reads are those where the number of complement events >= number of template events.

    # Column 1 = fast5_filename, 53 = template_called_events, 61 = complement_called_events, 100 = twod_fastq_seqlen
    callstats = np.loadtxt(incallstatspath, dtype='S100', usecols=(0, 52, 60, 99), delimiter='\t')
    hq = {}
    for i in range(1, len(callstats)):
        if int(callstats[i][3]) > 0:
            hq[callstats[i][0]] = int(callstats[i][1]) <= int(callstats[i][2])
        else:
            hq[callstats[i][0]] = 'ignore'

    # Iterate through each of the records in the fasta or fastq file and write the output to files
    # outdir/outprefix_2Dfilt.fasta/q and outdir/outprefix_2Drest.fasta/q

    outfilt_path = os.path.join(outdir, outprefix + '_2Dfilt.' + infiletype)
    outrest_path = os.path.join(outdir, outprefix + '_2Drest.' + infiletype)

    with open(outfilt_path, 'w') as filt_path:
        with open(outrest_path, 'w') as rest_path:
            for record in SeqIO.parse(in2Dallpath, infiletype):
                id = record.description
                filename = record.description.split('=')[-1]
                if hq[filename] == 'ignore':
                    continue
                if infiletype == 'fasta':
                    if hq[filename]:
                        file_path.write('>{0}\n{1}\n'.format(record.id, record.seq))
                    else:
                        rest_path.write('>{0}\n{1}\n'.format(record.id, record.seq))
                elif infiletype == 'fastq':
                    #seq = "\n".join([str(record.seq)[i:i+int(linewidth)] for i in range(0, len(record.seq), int(linewidth))])
                    bq = ''.join([chr(x+33) for x in record.letter_annotations['phred_quality']])
                    if hq[filename]:
                        filt_path.write('@{id}\n{seq}\n+\n{bq}\n'.format(id=id, seq=record.seq, bq=bq))
                    else:
                        rest_path.write('@{id}\n{seq}\n+\n{bq}\n'.format(id=id, seq=record.seq, bq=bq))

