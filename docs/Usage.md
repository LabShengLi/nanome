**This is an explanation of how to use NANOME pipeline on raw Fast5 input.**

The inputs of NANOME pipeline is a folder/tar/tar.gz or txt file list containing raw signal Fast5 files and a reference genome. We recommend allocate GPU resources to softwares such as Guppy, DeepSignal, DeepMod and Megalodon, in order to optimal running times. We integrated a latest tool METEORE's result, it is depends on other tools's read-level outputs (e.g., Megalodon and DeepSignal), and running METEORE program directly on them, detailed please check [METEORE](https://github.com/comprna/METEORE).

# 1. Running NANOME for human nanopore sequencing data


## Pre-defined pipeline profiles
`-profile` is the name of execution configuration, we support various of  configurations, e.g., `conda`, `docker`, `singularity`, `hpc` and `google`. Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Conda, Docker, Singularity) - see below.

Note that multiple profiles can be loaded, for example: `-profile singularity,hpc` - the order of arguments is important! They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the PATH. This is not recommended.

* `conda`
  * A generic configuration profile to be used with [Conda](https://docker.com/)
* `docker`
  * A generic configuration profile to be used with [Docker](https://docker.com/)
  * Pulls from [Docker Hub](https://hub.docker.com/repository/docker/liuyangzzu/nanome): liuyangzzu/nanome:latest
* `singularity`
  * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
  * Pulls from [Docker Hub](https://hub.docker.com/repository/docker/liuyangzzu/nanome): docker://liuyangzzu/nanome:latest
* `hpc`		
  * A generic configuration profile to be used on HPC cluster with [SLURM](https://slurm.schedmd.com/documentation.html) job submission support.
* `google`	
  * A generic configuration profile to be used on [Google Cloud](https://cloud.google.com/) platform with **google-lifesciences** support.

## Running samples
The command for running NANOME pipeline is to run `nextflow run https://github.com/TheJacksonLaboratory/nanome.git`. `--input` is input Fast5 file locations, our pipeline support three kinds of inputs: (1) folder, (2) tar/tar.gz file, (3) a txt file `.filelist.txt` contains list of compressed Fast5 files/folders. `--dsname` is output dataset name.

By default, we are using hg38 human reference genome, and you can specify other reference genome using parameter `dataType='ecoli'`. An example of how to use NANOME pipeline is given below.

```angular2html
# Get nextflow executable file
curl -fsSL get.nextflow.io | bash

# Get command help
nextflow run https://github.com/TheJacksonLaboratory/nanome.git --help

# Running NANOME pipeline on HPC using Singularity
nextflow run https://github.com/TheJacksonLaboratory/nanome.git \
    -profile ci,singularity,hpc

# Running NANOME pipeline for human data
nextflow run https://github.com/TheJacksonLaboratory/nanome.git \
    -profile singularity,hpc \
    --dsname TestData \
    --input https://raw.githubusercontent.com/TheJacksonLaboratory/nanome/master/inputs/test.demo.filelist.txt
```

You can also running NANOME pipeline on cloud computing platform ([google cloud platform](https://cloud.google.com/) or [Lifebit CloudOS](https://lifebit.gitbook.io/cloudos/)), sample of command line is below.
```angular2html
# Running on Google Cloud (https://cloud.google.com)
nextflow run https://github.com/TheJacksonLaboratory/nanome.git \
    -profile ci,docker,google \
    -w [Google-storage-bucket]/TestData-work \
    --outputDir [Google-storage-bucket]/TestData-ouputs
```

Pipeline running results is below, output directory trees are [outputs](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/outputs_demo.tree.txt) and [work](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/work_demo.tree.txt). It can also generates [timeline](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/timeline_demo.pdf), [report](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/report_demo.pdf) and [resource usage](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/trace_demo.txt.tsv) with more Nextflow [options](https://www.nextflow.io/docs/latest/tracing.html) (e.g., `-with-report -with-timeline -with-trace -with-dag -resume`).

```angular2html
N E X T F L O W  ~  version 20.10.0
Launching `main.nf` [wise_crick] - revision: efbaa90697
NANOME - NF PIPELINE (v1.3)
by Li Lab at The Jackson Laboratory
https://nanome.jax.org
=================================
dsname          :TestData
input           :https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt
output          :/fastscratch/li-lab/nanome/outputs
work            :/fastscratch/li-lab/nanome/work
dataType        :human
runBasecall     :true
runMethcall     :true
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


All tools's methlation calling and evaluation results will be output to `outputs` folder by default below.

```angular2html
tree outputs/TestData-methylation-callings/

outputs/TestData-methylation-callings/
├── Raw_Results-TestData
│   ├── TestData.deepmod.C_clusterCpG_per_site.combine.bed.gz
│   ├── TestData.deepmod.C_per_site.combine.bed.gz
│   ├── TestData.deepsignal.per_read.combine.tsv.gz
│   ├── TestData.guppy.fast5mod_per_site.combine.tsv.gz
│   ├── TestData.guppy.gcf52ref_per_read.combine.tsv.gz
│   ├── TestData.megalodon.per_read.combine.bed.gz
│   ├── TestData.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.tsv.gz
│   ├── TestData.nanopolish.per_read.combine.tsv.gz
│   └── TestData.tombo.per_read.combine.bed.gz
├── Read_Level-TestData
│   ├── TestData_DeepSignal-perRead-score.tsv.gz
│   ├── TestData_Guppy-perRead-score.tsv.gz
│   ├── TestData_Megalodon-perRead-score.tsv.gz
│   ├── TestData_METEORE-perRead-score.tsv.gz
│   ├── TestData_Nanopolish-perRead-score.tsv.gz
│   └── TestData_Tombo-perRead-score.tsv.gz
└── Site_Level-TestData
    ├── TestData_DeepMod-perSite-cov1.sort.bed.gz
    ├── TestData_DeepSignal-perSite-cov1.sort.bed.gz
    ├── TestData_Guppy-perSite-cov1.sort.bed.gz
    ├── TestData_Megalodon-perSite-cov1.sort.bed.gz
    ├── TestData_METEORE-perSite-cov1.sort.bed.gz
    ├── TestData_Nanopolish-perSite-cov1.sort.bed.gz
    └── TestData_Tombo-perSite-cov1.sort.bed.gz

tree outputs -L 1
outputs
├── README.txt
├── report
├── TestData-basecallings
└── TestData-methylation-callings
```

We also support input as a file list if input is suffix like `.filelist.txt`, an example input is [inputs/test.demo.filelist.txt](https://github.com/TheJacksonLaboratory/nanome/blob/master/inputs/test.demo.filelist.txt). Please use folowings for command line help:
```angular2html
nextflow run main.nf --help

N E X T F L O W  ~  version 20.10.0
Launching `main.nf` [maniac_babbage] - revision: deabf939f5
Usage:
The typical command for running the pipeline is as follows:
nextflow run main.nf -profile singularity,hpc --dsname TestData --input https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt

Mandatory arguments:
  --dsname                      Dataset name
  --input                       Input raw fast5 folders/tar/tar.gz files

Options:
  --processors                  Processors used for each process, default is 8
  --outputDir                   Output dir, default is 'outputs'
  --dataType                    Data type, default is 'human', can be also 'ecoli'
  --referenceGenome             Reference genome, default is 'reference_genome/hg38/hg38.fasta'

  --cleanCache                  True if clean work dir after complete
  --computeName                 Command used for tools, default is 'gpu', can be 'cpu'

  --queueName                   SLURM job submission queue name, default is 'gpu'
  --qosName                     SLURM job submission qos name, default is 'inference'
  --gresGPUOptions              SLURM job submission GPU options, default is '--gres=gpu:v100:1'
  --jobMaxTime              	SLURM job submission time options, default is '06:00:00'
  --jobMaxMem                   SLURM job submission memory options, default is '180G'

  --conda_name              	Conda name used for pipeline
  --docker_name              	Docker name used for pipeline
  --singularity_name            Singularity name used for pipeline
  --singularity_cache_dir       Singularity cache dir

Other options:
  --guppyDir                    Guppy installation dir

-profile options:
  Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

  docker 			A generic configuration profile to be used with Docker, pulls software from Docker Hub: quay.io/liuyangzzu/nanome:v1.5
  singulairy			A generic configuration profile to be used with Singularity, pulls software from: docker://quay.io/liuyangzzu/nanome:v1.5
  conda				Please only use conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity. Create conda enviroment by 'conda env create -f environment.yml'
  hpc				A generic configuration profile to be used on HPC with SLURM job submission support.
  google			A generic configuration profile to be used on Google Cloud with 'google-lifesciences' support.
```


# 2. Experiment for E. coli data
The NANOME pipeline supports 5mC detection by all tools on both human and Escherichia coli data. Note that `--dataType` need to be set as `ecoli`. Below is an example of pipeline runing on E. coli data, please refer to the input parameters for pipeline `-config` params [conf/ecoli_demo.config](https://github.com/TheJacksonLaboratory/nanome/blob/master/conf/ecoli_demo.config).

```angular2html
git clone https://github.com/TheJacksonLaboratory/nanome.git
cd nanome

nextflow run main.nf \
    -profile singularity,hpc \
    -config conf/ecoli_demo.config
```

Pipeline results for E. coli data is below.

```angular2html
N E X T F L O W  ~  version 20.10.0
Launching `main.nf` [drunk_carlsson] - revision: efbaa90697
NANOME - NF PIPELINE (v1.3)
by Li Lab at The Jackson Laboratory
https://nanome.jax.org
=================================
dsname          :EcoliDemo
input           :https://zenodo.org/record/5483859/files/ecoli_data_from_meteore.tar.gz
output          :/fastscratch/li-lab/nanome/outputs-ecoli
work            :/fastscratch/li-lab/nanome/work-ecoli
dataType        :ecoli
runBasecall     :true
runMethcall     :true
=================================
executor >  slurm (19)
[e6/9fbf7e] process > EnvCheck (EnvCheck)            [100%] 1 of 1 ✔
[10/24fad4] process > Untar (ecoli_data_from_mete... [100%] 1 of 1 ✔
[45/94553b] process > Basecall (ecoli_data_from_m... [100%] 1 of 1 ✔
[6d/7fd4ca] process > QCExport                       [100%] 1 of 1 ✔
[2c/31de5c] process > Resquiggle (ecoli_data_from... [100%] 1 of 1 ✔
[58/256829] process > Nanopolish (ecoli_data_from... [100%] 1 of 1 ✔
[67/bf7f23] process > Megalodon (ecoli_data_from_... [100%] 1 of 1 ✔
[97/dd7a8f] process > DeepSignal (ecoli_data_from... [100%] 1 of 1 ✔
[03/efbbca] process > Guppy (ecoli_data_from_mete... [100%] 1 of 1 ✔
[b9/3ed251] process > Tombo (ecoli_data_from_mete... [100%] 1 of 1 ✔
[4f/c40418] process > DeepMod (ecoli_data_from_me... [100%] 1 of 1 ✔
[5c/8182bf] process > NplshComb (1)                  [100%] 1 of 1 ✔
[d0/123535] process > MgldnComb (1)                  [100%] 1 of 1 ✔
[a1/25ef0c] process > DpSigComb (1)                  [100%] 1 of 1 ✔
[ff/f8786e] process > GuppyComb (1)                  [100%] 1 of 1 ✔
[07/e1fa8a] process > TomboComb (1)                  [100%] 1 of 1 ✔
[fa/863b04] process > DpmodComb (1)                  [100%] 1 of 1 ✔
[a0/3326d1] process > METEORE (1)                    [100%] 1 of 1 ✔
[2f/9a2ad7] process > Report (1)                     [100%] 1 of 1 ✔
Completed at: 11-Sep-2021 20:25:14
Duration    : 5m 19s
CPU hours   : 0.1
Succeeded   : 19
```

The output files of pipeline on E. coli data by all tools are below, please also check the pipeline output directory tree for [outputs](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/outputs_ecoli.tree.txt) and [work](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/work_ecoli.tree.txt). The pipeline can also generate [timeline](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/timeline_ecoli.pdf), [report](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/report_ecoli.pdf) and [resource usage](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/trace_ecoli.txt.tsv).


# 3. Benchmarking experiment
We constructed a list of benchmarking datasets that contain Fast5 reads from 800 to 8,000  for NA19240. The datasets can be got by users upon request. Following command is running NANOME pipeline on our benchmarking datasets, please refer to the input parameters for pipeline `-config` params [conf/benchmarking.config](https://github.com/TheJacksonLaboratory/nanome/blob/master/conf/benchmarking.config).

```angular2html
git clone https://github.com/TheJacksonLaboratory/nanome.git
cd nanome

nextflow run main.nf \
    -profile singularity,hpc  \
    -config conf/benchmarking.config
```

Resource usage are reported by [Nextflow](https://www.nextflow.io/) workflow reporting utilities. Please refer to the [Trace file](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/trace_benchmark.txt.tsv), [Report](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/report_benchmark.pdf) and [Timeline](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/timeline_benchmark.pdf) of benchmarking results on our HPC.

# 4. Running pipeline on cloud computing platform

Our Nextflow pipeline can running on CloudOS. The CloudOS will use a default Docker image. Below is an example.

```angular2html
git clone https://github.com/TheJacksonLaboratory/nanome.git
cd nanome

nextflow run main.nf \
    -profile ci,docker,google \
    -w [Google-storage-bucket]/nanome-work-ci \
    --outputDir [Google-storage-bucket]/nanome-outputs-ci
```

The `jax-nanopore-01-project-data` is a sample of **Data Bucket** name that you can access on google cloud. `-w` is pipeline output working directory, `--outputDir` is the directory for methylation-calling results.

For more detail of using cloud computing, please check [Cloud computing usage](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/CloudComputing.md).

We will update more examples here within a short time.