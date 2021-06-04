#!/bin/bash

set -x
date;hostname;pwd


sbatch --job-name=NA12878.chr1-3_f3.nxf.gcloud run_on_google_cloud.sbatch CHR1_3_f3 inputs/na12878.google_cloud.chr1-3_f3.filelist.txt

# sbatch --job-name=NA12878.chr1-3.nxf.gcloud run_on_google_cloud.sbatch CHR1_3 inputs/na12878.google_cloud.chr1-3.filelist.txt


#sbatch --job-name=NA12878.chr4-6.nxf.gcloud run_on_google_cloud.sbatch CHR4_6 inputs/na12878.google_cloud.chr4-6.filelist.txt

# sbatch --job-name=NA12878.chr20.nxf.gcloud run_on_google_cloud.sbatch CHR20 inputs/na12878.google_cloud.chr20.filelist.txt
