conda activate flye
date
flye --genome-size 1.1g --pacbio-raw pacbio.subreads.fa.gz --out-dir ./outdir/ --threads 30 --iterations 2 --resume
date

