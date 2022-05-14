# NANOME pipeline (Nanopore long-read sequencing data consensus DNA methylation detection)   

[![demo_gif.gif](https://github.com/LabShengLi/nanome/blob/master/docs/demo_gif.gif)](https://www.youtube.com/watch?v=TfotM55KTVE)

## Highlights of NANOME pipeline
### Several first highlights for NANOME

![Figure_pipe_comp](https://github.com/LabShengLi/nanome/blob/xgboost1/docs/resources/pipeline_comparison.jpg)

* Enables users to process **terabasescale** Oxford Nanopore sequencing datasets.
* Provide a **one command line**/**web-based UI** for end-to-end analyzing Nanopore sequencing methylation-callings.
* Support **various platform** executions: local, HPC and CloudOS, **without needs for tools' installation** (NANOME support docker and singularity).
* **First standardized whole genome-wide evaluation framework**, considering per-read and per-site performance for singletons/non-singletons, genic and intergenic regions, CpG islands/shores/shelves, different CG densities regions and repetitive regions. 
* The **first Nextflow based DNA methylation-calling pipeline for ONT data**. Please check more articles about Nextflow based workflow technology from Nature Biotechnology: https://doi.org/10.1038/s41587-020-0439-x and https://doi.org/10.1038/nbt.3820.
* Allow **add new modules/tools** in simple config txt file, without need to touch the main pipeline codes, supporting rapid development and evaluation.
* Consensus of top performers by XGBoost model, allow NA values.
* Multi-modifications for 5mC and 5hmC.
* Haplotype-awared phasing and allele-specific methylation detection.


## Methodology of NANOME pipeline

[comment]: <> (**Background:** Nanopore long-read sequencing technology greatly expands the capacity of long-range, single-molecule DNA-modification detection. A growing number of analytical tools have been developed to detect DNA methylation from nanopore sequencing reads. Here, we assess the performance of different methylation calling tools to provide a systematic evaluation to guide researchers performing human epigenome-wide studies.)


![Figure1A](https://github.com/LabShengLi/nanome/blob/master/docs/Fig1A.jpg)

**Fig. 1A. Survey of methylation calling tools .**  Timeline of publication and technological developments of Oxford Nanopore Technologies (ONT) methylation calling tools to detect DNA cytosine modifications. 


![Figure1B](https://github.com/LabShengLi/nanome/blob/master/docs/Fig1B.jpg)
**Fig. 1B. Workflow for 5-methylcytosine (5mC) detection for nanopore sequencing.** 


[comment]: <> (**Results:** We compared several analytic tools for detecting DNA modifications from nanopore long-read sequencing data. We evaluated the CpG methylation-detection accuracy, CpG site coverage, and running time using nanopore sequencing data across different genomic contexts, using natural human DNA. Furthermore, we provide an online DNA methylation database &#40;https://nanome.jax.org&#41; with which to display the DNA methylation levels detected by nanopore sequencing and bisulfite sequencing data across different genomic contexts.)


[comment]: <> (**Conclusions:** Our study is the first benchmark of state-of-the-art methods for detection of mammalian whole-genome DNA-modifications in nanopore sequencing. We provide a broad foundation for cross-platform standardization, and an evaluation of analytical tools designed for genome-scale modified-base detection using nanopore sequencing. )



### CI/CD features
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
Please refer to [Usage](https://github.com/LabShengLi/nanome/blob/master/docs/Usage.md) and [Specific Usage](https://github.com/LabShengLi/nanome/blob/master/docs/SpecificUsage.md) and [NANOME options](https://github.com/LabShengLi/nanome/blob/tutorial1/docs/nanome_params.md) for how to use NANOME pipeline. For running on CloudOS platform (e.g., google cloud), please check [Usage on CloudOS](https://github.com/LabShengLi/nanome/blob/master/docs/Usage.md#4-running-pipeline-on-cloud-computing-platform). We provide a **tutorial video** for running NANOME pipeline:

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
Please check [phasing usage](https://github.com/LabShengLi/nanome/blob/tutorial1/docs/Phasing.md).
![PhasingDemo](https://github.com/LabShengLi/nanome/blob/tutorial1/docs/resources/nanome3t_5mc_phasing2.png)

### Lifebit CloudOS report
We now support running NANOME on cloud computing platform. [Lifebit](https://lifebit.ai/lifebit-cloudos/) is a web-based cloud computing platform, and below is the running reports:
* Ecoli test report: https://cloudos.lifebit.ai/public/jobs/61c9fd328c574a01e8d31d2e
* Human test report: https://cloudos.lifebit.ai/public/jobs/61c9fe618c574a01e8d31e99
* NA12878 chr22 report: https://cloudos.lifebit.ai/public/jobs/61c4f2ad8c574a01e8d0eee3
* NA12878 chr20 part5 report: https://cloudos.lifebit.ai/public/jobs/61c770748c574a01e8d2062b


## Revision History
For release history, please visit [here](https://github.com/LabShengLi/nanome/releases). For details, please go [here](https://github.com/LabShengLi/nanome/blob/master/README.md).


## Contact
If you have any questions/issues/bugs, please post them on [GitHub](https://github.com/LabShengLi/nanome/issues). We will continuously update the GitHub to support famous methylation-calling tools for Oxford Nanopore sequencing.


## Reference

**DNA methylation-calling tools for Oxford Nanopore sequencing: a survey and human epigenome-wide evaluation.** Genome Biology 22, 295 (2021). https://doi.org/10.1186/s13059-021-02510-z  
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02510-z 
