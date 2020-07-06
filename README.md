Our comparative genomics work on covid-19.

Integrated into our dashboard at https://bioinfo.lau.edu.lb/gkhazen/covid19/

to complement the testing statistics timeseries (ref repo covid19-testing-data)


## How to use

The pipeline can be run as follows:

First set of files to process sequences from GISAID. This ultimately yields a VCF files containing the mutations:

- l1a_parse_gisaid_fasta.R, l1b_mummer.sh, l1c_stats.R (in this order)

The scripts assume the following files are in the same directory as the code:

- `Reference/EPI_ISL_402125.fasta`
- `MainFile/somefile.fasta` (that’s the file i download from gisaid)
- `orfs_covid19.txt`


Second set of files to calculate annotate with snpEff and prepare for upload to the dashboard

- l2a_snpEff.sh
- l2b_files4dashboard.R


Third script which postprocesses the dashboard file a bit more:

- `l3b_postprocess_genomics.py`: by Halim, used before uploading to arcgis.com
  - The jupyter notebook from which this was derived is the l3a file


