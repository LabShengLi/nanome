# NANOME pipeline (Nanopore sequencing consensus DNA methylation detection method and pipeline)   

[![demo_gif.gif](https://github.com/LabShengLi/nanome/blob/master/docs/demo_gif.gif)](https://www.youtube.com/watch?v=TfotM55KTVE)

## Highlights of NANOME pipeline
### Several first highlights for NANOME

![Figure_pipe_comp](https://github.com/LabShengLi/nanome/blob/master/docs/resources/pipeline_comparison.jpg)

* Enables users to process **terabasescale** Oxford Nanopore sequencing datasets.
* Provide a **one command line**/**web-based UI** for end-to-end analyzing Nanopore sequencing methylation-callings.
* Support **various platform** executions: local, HPC and CloudOS, **without needs for tools' installation** (NANOME support docker and singularity).
* **First standardized whole genome-wide evaluation framework**, considering per-read and per-site performance for singletons/non-singletons, genic and intergenic regions, CpG islands/shores/shelves, different CG densities regions and repetitive regions. 
* The **first Nextflow based DNA methylation-calling pipeline for ONT data**. Please check more articles about Nextflow based workflow technology from Nature Biotechnology: https://doi.org/10.1038/s41587-020-0439-x and https://doi.org/10.1038/nbt.3820.
* Allow **add new modules/tools** in simple config txt file, without need to touch the main pipeline codes, supporting rapid development and evaluation.
* Consensus of top performers by XGBoost model, allow NA values.
* Multi-modifications for 5mC and 5hmC.
* Haplotype-awared phasing and allele-specific methylation detection.
* Support [Dorado](https://github.com/nanoporetech/dorado) basecall and methylation call.


## Background

[comment]: <> (**Background:** Nanopore long-read sequencing technology greatly expands the capacity of long-range, single-molecule DNA-modification detection. A growing number of analytical tools have been developed to detect DNA methylation from nanopore sequencing reads. Here, we assess the performance of different methylation calling tools to provide a systematic evaluation to guide researchers performing human epigenome-wide studies.)


![Figure1A](https://github.com/LabShengLi/nanome/blob/master/docs/Fig1A.jpg)

**Survey of methylation calling tools .**  Timeline of publication and technological developments of Oxford Nanopore Technologies (ONT) methylation calling tools to detect DNA cytosine modifications. 


![Figure1B](https://github.com/LabShengLi/nanome/blob/master/docs/Fig1B.jpg)

**Workflow for 5-methylcytosine (5mC) detection for nanopore sequencing.** 


[comment]: <> (**Results:** We compared several analytic tools for detecting DNA modifications from nanopore long-read sequencing data. We evaluated the CpG methylation-detection accuracy, CpG site coverage, and running time using nanopore sequencing data across different genomic contexts, using natural human DNA. Furthermore, we provide an online DNA methylation database &#40;https://nanome.jax.org&#41; with which to display the DNA methylation levels detected by nanopore sequencing and bisulfite sequencing data across different genomic contexts.)


[comment]: <> (**Conclusions:** Our study is the first benchmark of state-of-the-art methods for detection of mammalian whole-genome DNA-modifications in nanopore sequencing. We provide a broad foundation for cross-platform standardization, and an evaluation of analytical tools designed for genome-scale modified-base detection using nanopore sequencing. )

## CI/CD automation features
We use  CI Automation Tools to **enable the automated testing on every commit and on PRs** to make sure that updates are not introducing bugs. Please check the automatic testing results on [Github](https://github.com/LabShengLi/nanome/actions).


## System Requirements
### Hardware requirements
NANOME pipeline can be easily configured with different RAM, CPU/GPU resources schema to parallelly run methylation-calling tools. For optimal usage, we recommend running NANOME pipeline on HPC or cloud computing platform, e.g., google cloud platform (GCP). The basic hardware requirements are below:
* GPU or CPU with 2+ cores. 
* RAM: 7+ GB per cpu.
* Storage using HDD or SSD. Please ensure the storage before running the pipeline.


### Software requirements
NANOME pipeline uses Nextflow technology. Users only need to install [Nextflow](https://www.nextflow.io/) (check the installation guide from https://nf-co.re/usage/installation), and have one of below commonly used environment tool:
* [Conda](https://docs.conda.io/en/latest/miniconda.html)
* [Docker](https://docs.docker.com/get-docker)
* [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html)

We provide conda, docker and singularity environments that depend on below well-known open-source packages for basecalling/methylation-calling/phasing on nanopore sequencing data:

[nanopolish](https://github.com/jts/nanopolish) >=0.13.2  
[megalodon](https://github.com/nanoporetech/megalodon) >=2.2.9  
[deepsignal](https://github.com/bioinfomaticsCSU/deepsignal) >=0.1.8  
[ont-tombo](https://github.com/nanoporetech/tombo) >=1.5.1  
[deepmod](https://github.com/WGLab/DeepMod) >=0.1.3  
[METEORE](https://github.com/comprna/METEORE) >=1.0.0  
[ont-pyguppy-client-lib](https://github.com/nanoporetech/pyguppyclient) >=4.2.2  
[fast5mod](https://github.com/nanoporetech/fast5mod) >=1.0.5  
[Clair3](https://github.com/HKU-BAL/Clair3) >=v0.1-r11  
[Whatshap](https://github.com/whatshap/whatshap) >=1.0  
[NanomethPhase bam2bis](https://github.com/vahidAK/NanoMethPhase) >= 1.0  
[GNU Parallel](https://www.gnu.org/software/parallel) >=20170422  


Guppy software >= 4.2.2 from [ONT (Oxford Nanopore Technologies) website](https://nanoporetech.com)


## Installation
Users only need to install **Nextflow** (https://nf-co.re/usage/installation). NANOME execution environment will be automatically configured with the support of conda, docker or singularity containers. Below is steps for installing Nextflow:
```angular2html
# Install nextflow
conda install -c conda-forge -c bioconda nextflow
nextflow -v
```

NANOME pipeline support running with various ways in different platforms:
* Docker
* Singularity
* Conda
* **Local** execution: running directly on default platform
* HPC clusters with **SLURM** support
* Cloud computing platform, e.g., Google Cloud Platform(GCP) with **google-lifesciences** support


## Simple usage
Please refer to [Usage](https://github.com/LabShengLi/nanome/blob/master/docs/Usage.md) and [Specific Usage](https://github.com/LabShengLi/nanome/blob/master/docs/SpecificUsage.md) and [NANOME options](https://github.com/LabShengLi/nanome/blob/master/docs/nanome_params.md) for how to use NANOME pipeline. For running on CloudOS platform (e.g., google cloud), please check [Usage on CloudOS](https://github.com/LabShengLi/nanome/blob/master/docs/Usage.md#5-running-pipeline-on-cloud-computing-platform). We provide a **tutorial video** for running NANOME pipeline:

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/TfotM55KTVE/0.jpg)](https://www.youtube.com/watch?v=TfotM55KTVE)

When you have Nextflow software, NANOME pipeline can be directly executed without any other additional installation steps:
```angular2html
# Run NANOME via docker
nextflow run LabShengLi/nanome\
    -profile test,docker

# Run NANOME via singularity
nextflow run LabShengLi/nanome\
    -profile test,singularity

# Run NANOME for human data
nextflow run LabShengLi/nanome\
    -profile test_human,[docker/singularity]
```
Please note that above commands are integrated in our **CI/CD test cases**. Our GitHub will automatically test and report results on every commit and PRs (https://github.com/LabShengLi/nanome/actions). 

We firstly proposed the **standardized whole genome-wide evaluation packages**, check [standardized evaluation tool usage](https://github.com/LabShengLi/nanome/blob/master/docs/Eval.md) for more detail. We do not suggest evaluating on a portion of CpGs for performance comparisons.


## Train and test script for consensus model in NANOME

We train an xgboost model on top performers: Nanopolish, DeepSignal and Megalodon, for detailed input/output format of consensus model train and predict, check [consensus model format](https://github.com/LabShengLi/nanome/blob/master/docs/XGBoost_format.md). 

The training script usage is below:
```angular2html
cs_train.py  -h
usage: cs_train (NANOME) [-h] --train TRAIN [TRAIN ...] --train-chr TRAIN_CHR
                         [TRAIN_CHR ...] --test TEST [TEST ...] --test-chr
                         TEST_CHR [TEST_CHR ...]
                         [--input-tools INPUT_TOOLS [INPUT_TOOLS ...]]
                         [--dsname DSNAME] [--model-name MODEL_NAME]
                         [--base-model BASE_MODEL] -o O [--niter NITER]
                         [--cv CV] [--scoring SCORING]
                         [--random-state RANDOM_STATE]
                         [--processors PROCESSORS] [--test-lines TEST_LINES]
                         [--show-confusion-matrix] [--apply-cutoff]
                         [--apply-cutoff-train] [--verbose]

Consensus model train on data

optional arguments:
  -h, --help            show this help message and exit
  --train TRAIN [TRAIN ...]
                        train data file
  --train-chr TRAIN_CHR [TRAIN_CHR ...]
                        train chr file
  --test TEST [TEST ...]
                        test data file
  --test-chr TEST_CHR [TEST_CHR ...]
                        train chr file
  --input-tools INPUT_TOOLS [INPUT_TOOLS ...]
                        input features for train, default is megalodon,
                        nanopolish, and deepsignal
  --dsname DSNAME       dataset name, default is NA12878
  --model-name MODEL_NAME
                        model name: basic, etc.
  --base-model BASE_MODEL
                        base model name: rf, xgboost, etc.
  -o O                  output file dir
  --niter NITER         number of iterations for random CV, default is 20
  --cv CV               number of CV, default is 3
  --scoring SCORING     optimized score name, i.e., f1, roc_auc, etc., default
                        is f1
  --random-state RANDOM_STATE
                        random state 42
  --processors PROCESSORS
                        number of processors, default is 1
  --test-lines TEST_LINES
                        test top N rows, such as 10000, default is None
  --show-confusion-matrix
                        if output verbose info
  --apply-cutoff        if apply default cutoff of tools
  --apply-cutoff-train  if apply default cutoff of tools before train
  --verbose             if output verbose info
```

Prediction script usage for xgboost model is below:
```angular2html
cs_predict.py -h
usage: cs_predict (NANOME) [-h] [-v] [-i I [I ...]] [--nanopolish NANOPOLISH]
                           [--megalodon MEGALODON] [--deepsignal DEEPSIGNAL]
                           [--feature FEATURE]
                           [--feature-readids-col FEATURE_READIDS_COL [FEATURE_READIDS_COL ...]]
                           [--feature-readids-col-order FEATURE_READIDS_COL_ORDER [FEATURE_READIDS_COL_ORDER ...]]
                           [--feature-seq-col FEATURE_SEQ_COL]
                           [--model_specific MODEL_SPECIFIC] -m M --dsname
                           DSNAME -o O [-t T [T ...]]
                           [--random-state RANDOM_STATE]
                           [--processors PROCESSORS] [--chunksize CHUNKSIZE]
                           [--inner-join] [--chrs CHRS [CHRS ...]]
                           [--interactive] [--verbose]

Consensus model predict for data

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -i I [I ...]          input tsv combined data for predicting
  --nanopolish NANOPOLISH
                        input nanopolish unified read-level file
  --megalodon MEGALODON
                        input megalodon unified read-level file
  --deepsignal DEEPSIGNAL
                        input deepsignal unified read-level file
  --feature FEATURE     input feature file for DNAseq
  --feature-readids-col FEATURE_READIDS_COL [FEATURE_READIDS_COL ...]
                        column index for ID, Chr, Pos and Strand
  --feature-readids-col-order FEATURE_READIDS_COL_ORDER [FEATURE_READIDS_COL_ORDER ...]
                        column index order for ID, Chr, Pos and Strand
  --feature-seq-col FEATURE_SEQ_COL
                        column index for DNA seq feature
  --model_specific MODEL_SPECIFIC
                        specific model info
  -m M                  model file, existing model list: NANOME2T,NANOME3T,xgb
                        oost_basic,xgboost_basic_w,xgboost_basic_w_seq
  --dsname DSNAME       dataset name
  -o O                  output file name
  -t T [T ...]          tools used for prediction, default is None
  --random-state RANDOM_STATE
                        random state, default is 42
  --processors PROCESSORS
                        num of processors, default is 8
  --chunksize CHUNKSIZE
                        chunk size for load large data, default is 500000
  --inner-join          if inner join for merge data, default is outer join
  --chrs CHRS [CHRS ...]
                        chromosomes used
  --interactive         if output to console as interactive mode, quit use q/Q
  --verbose             if output verbose info
```

Script for read-level performance comparison (accuracy, F1-score, etc.) on joined predictions by all tools:
```
cs_eval_read.py  -h

                        [--model-name MODEL_NAME [MODEL_NAME ...]]
                        [--model-file MODEL_FILE [MODEL_FILE ...]] -o O
                        [--processors PROCESSORS] [--bs-cov BS_COV]
                        [--tool-cov TOOL_COV] [--eval-type EVAL_TYPE]
                        [--model-base-dir MODEL_BASE_DIR]
                        [--test-lines TEST_LINES] [--chunksize CHUNKSIZE]
                        [--force-llr2] [--verbose]

Consensus model train on data

optional arguments:
  -h, --help            show this help message and exit
  -i I [I ...]          input data file
  --dsname DSNAME       dataset name, default is NA12878
  --model-name MODEL_NAME [MODEL_NAME ...]
                        model name: rf, xgboost, etc.
  --model-file MODEL_FILE [MODEL_FILE ...]
                        model file
  -o O                  output file dir
  --processors PROCESSORS
                        number of processors, default is 1
  --bs-cov BS_COV       bs-seq coverage cutoff, default is 5
  --tool-cov TOOL_COV   ONT tool coverage cutoff, default is 1
  --eval-type EVAL_TYPE
                        evaluation type, read-level or site-level
  --model-base-dir MODEL_BASE_DIR
                        model file's base dir
  --test-lines TEST_LINES
                        test top N rows, such as 10000, default is None
  --chunksize CHUNKSIZE
                        chunk size for load large data, default is 500000
  --force-llr2          if convert megalodon llr to llr2
  --verbose             if output verbose info
```

Script for site-level performance comparison (MSE, PCC) on joined predictions by all tools:
```
cs_eval_site.py -h

                        [--processors PROCESSORS] [--bs-cov BS_COV]
                        [--tool-cov TOOL_COV] [--eval-type EVAL_TYPE]
                        [--model-base-dir MODEL_BASE_DIR]
                        [--test-lines TEST_LINES] [--chunksize CHUNKSIZE]
                        [--save-data SAVE_DATA] [--force-llr2] [--verbose]

Consensus model train on data

optional arguments:
  -h, --help            show this help message and exit
  -i I [I ...]          input data file
  --dsname DSNAME       dataset name, default is NA12878
  --model-name MODEL_NAME [MODEL_NAME ...]
                        model name: rf, xgboost, etc.
  --model-file MODEL_FILE [MODEL_FILE ...]
                        model file
  -o O                  output file dir
  --processors PROCESSORS
                        number of processors, default is 1
  --bs-cov BS_COV       bs-seq coverage cutoff, default is 5
  --tool-cov TOOL_COV   ONT tool coverage cutoff, default is 1
  --eval-type EVAL_TYPE
                        evaluation type, i.e., site-level
  --model-base-dir MODEL_BASE_DIR
                        model file's base dir
  --test-lines TEST_LINES
                        test top N rows, such as 10000, default is None
  --chunksize CHUNKSIZE
                        chunk size for load large data, default is 500000
  --save-data SAVE_DATA
                        if save prediction outputs
  --force-llr2          if convert megalodon llr to llr2
  --verbose             if output verbose info
```


## Pipeline reports for NANOME
### Benchmarking reports on our HPC using [Nextflow](https://www.nextflow.io/)
We constructed a set of benchmarking datasets that contain reads from 800 to about 7,200 reads for NA19240, and monitored job running timeline and resource usage on our HPC, reports generated by **Nextflow** workflows are: [Trace file](https://github.com/LabShengLi/nanome/blob/master/docs/resources/trace_benchmark.txt.tsv), [Report](https://github.com/LabShengLi/nanome/blob/master/docs/resources/report_benchmark.pdf)  and [Timeline](https://github.com/LabShengLi/nanome/blob/master/docs/resources/timeline_benchmark.pdf). 

Our HPC hardware specifications are as follows:
* CPU: Intel(R) Xeon(R) Gold 6136 CPU @ 3.00GHz
* GPU: Tesla V100-SXM2-32GB 
* RAM: 300 GB
* Slurm manager version: 19.05.5

Timeline figure for benchmarking experiments are below:
![Bench-timeline](https://github.com/LabShengLi/nanome/blob/master/docs/resources/timeline_benchmark.jpg)


### Pipeline DAG
![NanomeDag](https://github.com/LabShengLi/nanome/blob/master/docs/nanome_dag.png)


### NANOME report
Please check [NANOME report](https://github.com/LabShengLi/nanome/blob/master/docs/NANOME_report_html.pdf) for the sample report by NANOME pipeline.

![NanomeReportHtml](https://github.com/LabShengLi/nanome/blob/master/docs/nanome_report_html.png)


### Haplotype-aware consensus methylations
Please check [phasing usage](https://github.com/LabShengLi/nanome/blob/master/docs/Phasing.md).
![PhasingDemo](https://github.com/LabShengLi/nanome/blob/master/docs/resources/nanome3t_5mc_phasing2.png)

### Lifebit CloudOS report
We now support running NANOME on cloud computing platform. [Lifebit](https://lifebit.ai/lifebit-cloudos/) is a web-based cloud computing platform, and below is the running reports:
* Ecoli test report: https://cloudos.lifebit.ai/public/jobs/6430509445941801546e5f8f
* Human test report: https://cloudos.lifebit.ai/public/jobs/6430639045941801546e627f
* NA12878 chr22 report: https://cloudos.lifebit.ai/public/jobs/6430b64645941801546e7400


## Revision History
For release history, please visit [here](https://github.com/LabShengLi/nanome/releases). For details, please go [here](https://github.com/LabShengLi/nanome/blob/master/README.md).


## Contact
If you have any questions/issues/bugs, please post them on [GitHub](https://github.com/LabShengLi/nanome/issues). We will continuously update the GitHub to support famous methylation-calling tools for Oxford Nanopore sequencing.


[//]: # (## Reference)

[//]: # ()
[//]: # (**DNA methylation-calling tools for Oxford Nanopore sequencing: a survey and human epigenome-wide evaluation.** Genome Biology 22, 295 &#40;2021&#41;. https://doi.org/10.1186/s13059-021-02510-z )

