**This is an explanation of how to use NANOME pipeline on raw Fast5 input. For specific scenarios, please check [Specific Usage](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/SpecificUsage.md).**

The inputs of NANOME pipeline is a folder/tar/tar.gz or txt file list containing raw signal Fast5 files and a reference genome. We recommend allocate GPU resources to softwares such as Guppy, DeepSignal, DeepMod and Megalodon, in order to optimal running times. We integrated a latest tool METEORE, it depends on other tools' read-level outputs (e.g., Megalodon and DeepSignal), and running METEORE program directly on them, detailed please check [METEORE](https://github.com/comprna/METEORE).

# 1. Running NANOME for human nanopore sequencing data

## Running samples
The command for running NANOME pipeline is to run `nextflow run TheJacksonLaboratory/nanome`.
- `--dsname` is dataset/analysis name.
- `--input` is input Fast5 files path. Nanome pipeline support three kinds of inputs: (1) folder, (2) tar/tar.gz file, (3) a txt file `.filelist.txt` contains list of compressed Fast5 files/folders.
- `--genome` is reference genome.

By default, we are using `--genome=hg38` for human reference genome, and you can specify other reference genome using parameter `--genome=ecoli`. We defined a bunch of predefined running configuration params in profile in next section. An example of how to use NANOME pipeline is given below.

```angular2html
# Get pipeline help
nextflow run TheJacksonLaboratory/nanome --help

# Running NANOME pipeline for human data on HPC
nextflow run TheJacksonLaboratory/nanome\
    -profile singularity,hpc\
    --dsname TestData\
    --input https://github.com/TheJacksonLaboratory/nanome/raw/master/test_data/demo1_fast5_reads.tar.gz\
    --genome hg38\
    --queue gpu --qos inference --memory 32GB --time 1h --gresOptions gpu:v100:1

# Running NANOME pipeline for E. coli data on HPC
nextflow run TheJacksonLaboratory/nanome\
    -profile singularity,hpc\
    --dsname EcoliData\
    --input https://storage.googleapis.com/jax-nanopore-01-project-data/nanome-input/ecoli_data_from_meteore.tar.gz\
    --genome ecoli\
    --queue gpu --qos inference --memory 32GB --time 1h --gresOptions gpu:v100:1
```

## Methylation-calling tool configuration
By default, NANOME pipeline will execute top four performers: **Nanopolish, Megalodon, DeepSignal and Guppy**, we also provide a **NANOME concensus results** using XGBoost model trained on all fully methylated and unmethylated CpGs based on Nanopolish, Megalodon and DeepSignal outputs. The model is very robust and can deal with NA values, and can make prediction if there is a prediction by any tool.

The NANOME concensus results can cover more CpGs than any single tool, and perform a slightly better performance. User can supply params `--run[tool-name]` with values `true` or `false` to configure if running a specific tool. 

## Pre-defined pipeline profiles
`-profile` is the name of execution configuration, we support various of  configurations, e.g., `conda`, `docker`, `singularity`, `hpc` and `google`. Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Conda, Docker, Singularity) - see below.

Note that multiple profiles can be loaded, for example: `-profile singularity,hpc` - the order of arguments is important! They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the PATH. This is not recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](https://docker.com/)
  * Pulls from [Docker Hub](https://hub.docker.com/repository/docker/liuyangzzu/nanome): liuyangzzu/nanome:latest
* `singularity`
  * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
  * Pulls from [Docker Hub](https://hub.docker.com/repository/docker/liuyangzzu/nanome): docker://liuyangzzu/nanome:latest
* `conda`
  * A generic configuration profile to be used with [Conda](https://docker.com/), check [conda usage](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/Usage.md#5-conda-environment-for-local-running)
* `hpc`		
  * A generic configuration profile to be used on HPC cluster with [SLURM](https://slurm.schedmd.com/documentation.html) job submission support.
* `google`	
  * A generic configuration profile to be used on [Google Cloud](https://cloud.google.com/) platform with **google-lifesciences** support.

You can also running NANOME pipeline on cloud computing platform ([google cloud platform](https://cloud.google.com/) or [Lifebit CloudOS](https://lifebit.gitbook.io/cloudos/)), sample of command line is below.
```angular2html
# Running test on Google Cloud (https://cloud.google.com)
nextflow run TheJacksonLaboratory/nanome\
    -profile test,docker,google \
    -w [Google-storage-bucket]/TestData-work \
    --outputDir [Google-storage-bucket]/TestData-ouputs\
    --googleProjectName  [Google-project-name]
```

## Running results and outputs

Pipeline running results is below, output directory trees are [outputs](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/outputs_demo.tree.txt) and [work](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/work_demo.tree.txt). It can also generates [timeline](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/timeline_demo.pdf), [report](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/report_demo.pdf) and [resource usage](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/trace_demo.txt.tsv) with more Nextflow [options](https://www.nextflow.io/docs/latest/tracing.html) (e.g., `-with-report -with-timeline -with-trace -with-dag -resume`).

```angular2html
N E X T F L O W  ~  version 20.10.0
Launching `main.nf` [wise_crick] - revision: efbaa90697
NANOME - NF PIPELINE (v1.3.6)
by Li Lab at The Jackson Laboratory
https://nanome.jax.org
=================================
executor >  slurm (28)
[30/56fb5a] process > EnvCheck (EnvCheck)            [100%] 1 of 1 ✔
[79/2d4393] process > Untar (demo1.fast5.reads.tar)  [100%] 2 of 2 ✔
[c2/1e6506] process > Basecall (demo1.fast5.reads... [100%] 2 of 2 ✔
[fc/3f5bd2] process > QCExport                       [100%] 1 of 1 ✔
[61/a78829] process > Resquiggle (demo1.fast5.rea... [100%] 2 of 2 ✔
[ef/ce1142] process > Nanopolish (demo1.fast5.rea... [100%] 2 of 2 ✔
[cb/72762f] process > Megalodon (demo1.fast5.read... [100%] 2 of 2 ✔
[2c/217ece] process > DeepSignal (demo1.fast5.rea... [100%] 2 of 2 ✔
[b7/7140d1] process > Guppy (demo1.fast5.reads.tar)  [100%] 2 of 2 ✔
[04/1b4c2a] process > Tombo (demo1.fast5.reads.tar)  [100%] 2 of 2 ✔
[d9/43d3cf] process > DeepMod (demo1.fast5.reads.... [100%] 2 of 2 ✔
[5a/244b39] process > NplshComb (1)                  [100%] 1 of 1 ✔
[72/a4f1a7] process > MgldnComb (1)                  [100%] 1 of 1 ✔
[a8/151a61] process > DpSigComb (1)                  [100%] 1 of 1 ✔
[d3/0d1141] process > GuppyComb (1)                  [100%] 1 of 1 ✔
[83/c7ecd3] process > TomboComb (1)                  [100%] 1 of 1 ✔
[63/5e39c1] process > DpmodComb (1)                  [100%] 1 of 1 ✔
[09/a1d741] process > METEORE (1)                    [100%] 1 of 1 ✔
[10/e307c2] process > Report (1)                     [100%] 1 of 1 ✔
Completed at: 11-Sep-2021 20:16:38
Duration    : 15m 24s
CPU hours   : 1.0
Succeeded   : 28
```


All tools' methlation calling and evaluation results will be output to `outputs` folder by default below.

```angular2html
tree results/TestData-methylation-callings/

results/TestData-methylation-callings/
├── Raw_Results-TestData
│   ├── TestData_deepmod_clusterCpG_per_site_combine.bed.gz
│   ├── TestData_deepmod_c_per_site_combine.bed.gz
│   ├── TestData_deepsignal_per_read_combine.tsv.gz
│   ├── TestData_guppy_fast5mod_per_site_combine.tsv.gz
│   ├── TestData_megalodon_per_read_combine.bed.gz
│   ├── TestData_meteore_deepsignal_megalodon_optimized_rf_model_per_read_combine.tsv.gz
│   ├── TestData_nanome_NA12878_XGBoostNA3T_per_read_combine.tsv.gz
│   ├── TestData_nanopolish_per_read_combine.tsv.gz
│   └── TestData_tombo_per_read_combine.bed.gz
├── Read_Level-TestData
│   ├── TestData_DeepSignal-perRead-score.tsv.gz
│   ├── TestData_Megalodon-perRead-score.tsv.gz
│   ├── TestData_METEORE-perRead-score.tsv.gz
│   ├── TestData_NANOME-perRead-score.tsv.gz
│   ├── TestData_Nanopolish-perRead-score.tsv.gz
│   └── TestData_Tombo-perRead-score.tsv.gz
├── Site_Level-TestData
│   ├── TestData_DeepMod-perSite-cov1.sort.bed.gz
│   ├── TestData_DeepSignal-perSite-cov1.sort.bed.gz
│   ├── TestData_Guppy-perSite-cov1.sort.bed.gz
│   ├── TestData_Megalodon-perSite-cov1.sort.bed.gz
│   ├── TestData_METEORE-perSite-cov1.sort.bed.gz
│   ├── TestData_NANOME-perSite-cov1.sort.bed.gz
│   ├── TestData_Nanopolish-perSite-cov1.sort.bed.gz
│   └── TestData_Tombo-perSite-cov1.sort.bed.gz
└── tools_version_table.tsv

tree results -L 1
results
├── TestData-basecallings
├── TestData-methylation-callings
└── TestData_nanome_report.html
```

We also support input as a file list if input file name is suffixed like `.filelist.txt`, an example input is [test.demo.filelist.txt](https://github.com/TheJacksonLaboratory/nanome/blob/master/inputs/test.demo.filelist.txt). Please use folowings for pipeline command help:
```angular2html
nextflow run TheJacksonLaboratory/nanome --help
```

# 2. Experiment for E. coli data
The NANOME pipeline supports 5mC detection by all tools on both human and Escherichia coli data. Note that `--genome` need to be set as `ecoli`. Below is an example of pipeline runing on E. coli data, please refer to the input parameters for pipeline params' config file [ecoli_demo.config](https://github.com/TheJacksonLaboratory/nanome/blob/master/conf/examples/ecoli_demo.config).

```angular2html
nextflow run TheJacksonLaboratory/nanome\
    -profile singularity,hpc \
    -config  conf/examples/ecoli_demo.config
```

Pipeline results for E. coli data is below.

```angular2html
N E X T F L O W  ~  version 20.10.0
Launching `main.nf` [maniac_poitras] - revision: 47f69be0ab
NANOME - NF PIPELINE (v1.3.6)
by Li Lab at The Jackson Laboratory
https://github.com/TheJacksonLaboratory/nanome
=================================
executor >  slurm (14)
[61/6ba4f5] process > EnvCheck (EnvCheck)                      [100%] 1 of 1 ✔
[82/f5adec] process > Untar (ecoli_data_from_meteore.tar)      [100%] 1 of 1 ✔
[7b/4b3001] process > Basecall (ecoli_data_from_meteore.tar)   [100%] 1 of 1 ✔
[e8/4be95e] process > QCExport (EcoliDemo)                     [100%] 1 of 1 ✔
[49/000f5c] process > Resquiggle (ecoli_data_from_meteore.tar) [100%] 1 of 1 ✔
[cf/149a24] process > Nanopolish (ecoli_data_from_meteore.tar) [100%] 1 of 1 ✔
[9c/edf13d] process > NplshComb (EcoliDemo)                    [100%] 1 of 1 ✔
[e3/f56451] process > Megalodon (ecoli_data_from_meteore.tar)  [100%] 1 of 1 ✔
[87/9417c5] process > MgldnComb (EcoliDemo)                    [100%] 1 of 1 ✔
[f6/5ffa6b] process > DeepSignal (ecoli_data_from_meteore.tar) [100%] 1 of 1 ✔
[3f/b28b4d] process > DpSigComb (EcoliDemo)                    [100%] 1 of 1 ✔
[b0/93b421] process > Guppy (ecoli_data_from_meteore.tar)      [100%] 1 of 1 ✔
[db/f8cce4] process > GuppyComb (EcoliDemo)                    [100%] 1 of 1 ✔
[1d/5aba6f] process > Report (EcoliDemo)                       [100%] 1 of 1 ✔
Completed at: 22-Dec-2021 11:39:24
Duration    : 4m 54s
CPU hours   : 0.7
Succeeded   : 14
```

The output files of pipeline on E. coli data by all tools are below, please also check the pipeline output directory tree for [outputs](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/outputs_ecoli.tree.txt) and [work](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/work_ecoli.tree.txt). The pipeline can also generate [timeline](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/timeline_ecoli.pdf), [report](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/report_ecoli.pdf) and [resource usage](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/trace_ecoli.txt.tsv).

# 3. Support for other reference genome
We now support other reference genome. Below is an example of running NANOME for any other reference genomes, please make sure you put reference genome file .fasta and the indexed file into directory [reference-genome-dir], the `--chrSet` is the chomosomes params for the specific genome. 

```angular2html
nextflow run TheJacksonLaboratory/nanome\
    -profile singularity \
    --dsname [your-dataset-name]\
    --input [input-file]\
    --genome [reference-genome-dir]\
    --chrSet '[chomosomes sperated by a space]'
```

Note: NANOME support the default behaviour of each tool's running, if the tool perfomed on the genome you specified is same as human/E. coli data. 

# 4. Benchmarking experiment
We constructed a list of benchmarking datasets that contain Fast5 reads from 800 to 7,200  for NA19240. The datasets can be got by users upon request. Following command is running NANOME pipeline on our benchmarking datasets, please refer to the input parameters for config file [benchmarking_hpc.config](https://github.com/TheJacksonLaboratory/nanome/blob/master/conf/executors/benchmarking_hpc.config).

```angular2html
nextflow run TheJacksonLaboratory/nanome\
    -profile singularity,hpc  \
    -config  conf/executors/benchmarking_hpc.config\
    --dsname  BenchmarkData\
    --input  inputs/benchmark.filelist.txt
```

Resource usage are reported by [Nextflow](https://www.nextflow.io/) workflow reporting utilities. Please refer to the [Trace file](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/trace_benchmark.txt.tsv), [Report](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/report_benchmark.pdf) and [Timeline](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/timeline_benchmark.pdf) of benchmarking results on our HPC.

# 5. Running pipeline on cloud computing platform

Our Nextflow pipeline can running on CloudOS. The CloudOS recommend using the Docker image. Below is an example.

```angular2html
nextflow run TheJacksonLaboratory/nanome\
    -profile test,docker,google \
    -w [Google-storage-bucket]/nanome-work-ci \
    --outdir [Google-storage-bucket]/nanome-outputs-ci\
    --googleProjectName  [Google-project-name]
```

The `[Google-project-name]` is your google project name, and `[Google-storage-bucket]` is the **Data Bucket** name that you can access on google cloud. `-w` is pipeline output working directory, `--outdir` is the directory for methylation-calling results.

For more detail of using cloud computing, please check [Cloud computing usage](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/CloudComputing.md).

# 6. Conda environment for local running

NANOME support local running without Docker or Singularity support. Below is conda environment installation steps, users need to install Guppy software by themselves in this case:
```angular2html
# Create conda environment for local running NANOME
git clone https://github.com/TheJacksonLaboratory/nanome.git
cd nanome
conda env create --name nanome --file=environment.yml
conda activate nanome

pip install megalodon==2.3.5
npm install -g inliner

# Run NANOME pipeline using local execution
conda activate nanome
nextflow run TheJacksonLaboratory/nanome\
    -profile test\
    --guppyDir [guppy-installation-directory]

# Run NANOME pipeline via conda environment
nextflow run TheJacksonLaboratory/nanome\
    -profile test,conda\
    --conda_name [conda-env-dir]\
    --conda_base_dir [conda-dir]\
    --guppyDir [guppy-installation-directory]
```
Param`--guppyDir=[guppy-installation-directory]` is the Guppy software installation base directory, `--conda_base_dir [conda-dir]` is conda software base directory, `--conda_name [conda-env-dir]` is conda environment base directory.
