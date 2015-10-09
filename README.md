#poreQC v0.2.10

poreQC nanopore sequencing experiment report

Produces a PDF report based on either the intermediate or final output
data from an Oxford Nanopore MinION experiment.

Input:
* a file of meta-data related to sample provenance, library preparation and
  suitable reference for the target genome(s) with which to do the QC.
* directory containing a copy of the basecalled and/or pre-basecalled .fast5 files
  produced by the run so far.

Output:
* a directory containing one text file of tabular data required for each plot
* a directory of plots with suitable size and resolution for inclusion in a PDF
* a printable PDF version of the report

Notes:
# The program was designed to be run after an rsync of the raw sequencing data
   from the sequencing laptop to a server which can process the data without
   affecting the computer controlling the sequencing experiment. A special
   'data status file' will indicate whether the run has completed or is still ongoing.

```
usage:

  poreqc.py [-h] --inrundir str [--metadatafile str] [--statusfile str]
                 --refdir str --outdir str --outprefix str [--bconly]
                 [--forcereadstep] [--forcecallstep] [--forcereportstep]
                 [--overwrite] [--readlenbucketsz int] [--imgszconst int]
                 [--debug]

optional arguments:
  -h, --help            show this help message and exit
  --inrundir str        Top-level directory containing a copy of the data from
                        a sequencing run (default: None)
  --metadatafile str    Filename of the metadata file that should be in
                        --outdir (default: metadata.txt)
  --statusfile str      Filename of the status file that should be in --outdir
                        (Must contain either "Running" or "Finished")
                        (default: status.txt)
  --refdir str          Top-level directory for all references (e.g., for
                        wtchgR00000021 you would specify
                        /bsgdata/microbial/data/references to find file
                        wtchgR00000021/wtchgR00000021.fasta) (default:
                        /bsgdata/microbial/data/references)
  --outdir str          Write all output to this directory (default: None)
  --outprefix str       All output files will start with this prefix (default:
                        None)
  --bconly              Extract all information from the basecalled fast5
                        files only (default: False)
  --forcereadstep       Force recompute of readstats.txt and associated
                        summary tables (default: False)
  --forcecallstep       Force recompute of callstats.txt (default: False)
  --forcereportstep     Force recompute of PDF report and associated images
                        and summary statistics (default: False)
  --overwrite           Force recomputation of all output files, overwriting
                        existing as necessary (default: False)
  --readlenbucketsz int
                        Size of read length buckets on x-axis of MinKNOW-like
                        read yield distribution plot (default: 1000)
  --imgszconst int      Size parameter for images (default: 100)
  --debug               Print verbose diagnostics (value must be 0 or 1)
                        (default: False)
```
