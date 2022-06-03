## Tutorial of DNA methylation calling for ONT data
In this tutorial (20 min ~ 30 min), you will learn how to perform methylation calling on Oxford Nanopore sequencing data by latest tools. Please create a freshing new folder to execute following commands.


**Prerequisites: Assume your system has (1) basic command line utils (`curl`, `wget`, `zcat` and `tar`); (2) container supported by `Singularity` or `Docker`; (3) good internet connection. (Ability to connect to JAX Sumner HPC will accelerate working through the tutorial.)**

**This tutorial is compatible with different platforms: Linux, Mac, Windows, and CloudOS.** 


> ## Notes:
> **If you are HPC users, enter into an interactive node before running pipeline.**

> Enter an interactive node with 8 cpus for parallelly job running (HPC users only):

> ```
> srun --pty -q batch --time=01:00:00 --mem=25G -n 8  bash
> ```

## 1. Software installations
### 1.1 Install Docker or Singularity
In this tutorial, the ONT developed basecalling tool [Guppy](https://community.nanoporetech.com), methylation-calling tool [Megalodon](https://github.com/nanoporetech/megalodon) and our [NANOME](https://github.com/LabShengLi/nanome) consensus methylation detection pipeline use containerized environment supported by Docker or Singularity. 

Container environment avoids users to encounter software installations steps/issues, and ensures running belowing same commands across different platforms (Linux, MacOS and Windows). **If your system already has Singularity or Docker container, please skip this section.**

Install Docker from here: [https://docs.docker.com/get-docker](https://docs.docker.com/get-docker).

Install Singularity from here: [https://sylabs.io/guides/3.0/user-guide/installation.html](https://sylabs.io/guides/3.0/user-guide/installation.html).

#### 1.2 Install Nextflow
Running [NANOME](https://github.com/LabShengLi/nanome) consensus methylation detection pipeline needs install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html).

* If you have [Java 11 or later](https://www.oracle.com/java/technologies/downloads), you can install Nextflow with below command:


```
curl -s https://get.nextflow.io | bash

nextflow -v
```

* If you do not have Java, you can firstly install Conda, please follow this link ([https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)) to install Conda, such as [Install on Linux](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html):

```
wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Then you can install Nextflow through Conda:

```
conda create --name nanome python=3.9
conda activate nanome

# Install nextflow
conda install -c conda-forge -c bioconda nextflow

nextflow -v
```

### 2. 5mC & 5hmC detection by ONT developed tool Megalodon
In this section, you will prepare the Oxford Nanopore raw FAST5 files and a reference genome for input data.

#### 2.1 Download FAST5 files

```
wget https://github.com/LabShengLi/nanome/raw/master/test_data/ecoli_ci_test_fast5.tar.gz &&\
    tar -xzf ecoli_ci_test_fast5.tar.gz
    
ls ecoli_ci_test_fast5/
```

#### 2.2 Download reference genome

```
wget https://storage.googleapis.com/jax-nanopore-01-project-data/nanome-input/ecoli.tar.gz &&\
    tar -xzf ecoli.tar.gz

ls ecoli/
```

#### 2.3 5mC & 5hmC detection by ONT developed tool Megalodon
* Define a bash variable for Docker running (**Docker user only**):


```
RUN_NANOME="docker run -v $PWD:$PWD -w $PWD -it liuyangzzu/nanome"
```


* Define a bash variable for Singularity running (**Singularity user only**):


```
RUN_NANOME="singularity exec -e docker://liuyangzzu/nanome"
```


> ## Notes:
> For HPC users, `Singularity` may need to be loaded from module:
> ```
> module load singularity
> singularity  --version
> ```

* Check Guppy basecalling tool and Megalodon methylation-calling tool versions in container (**Note: The first time execution will take times (15 mins ~ 30 mins) to cache containers from DockerHub. The next time execution will be very quick.**):


```
$RUN_NANOME guppy_basecaller -v

$RUN_NANOME megalodon -v
```

Everything looks good now. You can run below command to perform basecalls, mappings, and CpG 5mC and 5hmC methylation-calls in both per-read (``mod_mappings``) and aggregated (``mods``) formats on prepared example ONT data.

```
$RUN_NANOME   megalodon \
    ecoli_ci_test_fast5/ \
    --guppy-config  dna_r9.4.1_450bps_fast.cfg\
    --remora-modified-bases dna_r9.4.1_e8 fast 0.0.0 5hmc_5mc CG 0 \
    --outputs basecalls mappings mod_mappings mods per_read_mods \
    --guppy-server-path  guppy_basecall_server \
    --reference ./ecoli/Ecoli_k12_mg1655.fasta \
    --output-directory methcall_ecoli_data\
    --write-mods-text --overwrite 

ls methcall_ecoli_data/
```

> ## Notes:
> **If your system supports GPU, the option `--devices 0` can be used for acceleration.**

> For meanings of options in Megalodon, please check link [https://github.com/nanoporetech/megalodon#getting-started](https://github.com/nanoporetech/megalodon#getting-started).

### 3. Consensus 5mC detection by NANOME Nextflow pipeline

#### 3.1 Run NANOME pipeline

We developed NANOME, the first Nextflow based pipeline for consensus DNA methylation detection using XGBoost, a gradient boosting algorithm for nanopore long-read sequencing. The consensus outputs can obtain more accurate performance (9%-13% MSE improvement) and  comprehensive CpG coverage (1%-7% more CpGs). NANOME pipeline supports input from both local locations and internet/cloud storages.

For execution of NANOME consensus pipeline for methylation detection, if you use Singularity container, specify `-profile singularity`; for Docker container, use `-profile docker` instead. Below is an example of using Singularity container on JAX Sumner HPC:

```
nextflow run LabShengLi/nanome\
    -profile singularity\
    --dsname CIEcoli\
    --input  https://github.com/LabShengLi/nanome/raw/master/test_data/ecoli_ci_test_fast5.tar.gz \
    --genome https://storage.googleapis.com/jax-nanopore-01-project-data/nanome-input/ecoli.tar.gz
```


`-profile` is the bundle of paramters, `--dsname` is the dataset name, `--input` is the FAST5 input files, and `--genome` is the reference genome file. For more options of NANOME, please check [NANOME options](https://github.com/LabShengLi/nanome/blob/tutorial1/docs/nanome_params.md)


> ## Notes:
> For Mac/PC users with low bandwidths/firewall permissions for accessing internet network, below is the local input running script:
> ```
> # Download DeepSignal model input file
> wget https://storage.googleapis.com/jax-nanopore-01-project-data/nanome-input/model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+.tar.gz
> 
> # Run NANOME on local input file
> nextflow run LabShengLi/nanome \
>   -profile docker \
>   --dsname CIEcoli \
>   --input ecoli_ci_test_fast5.tar.gz \
>   --genome ecoli.tar.gz \
>   --deepsignalDir model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+.tar.gz
> ```


The output of NANOME pipeline is followings:

```
N E X T F L O W  ~  version 22.04.0
Launching `https://github.com/LabShengLi/nanome` [stupefied_keller] DSL2 - revision: bcbf29ed93 [master]
NANOME - NF PIPELINE (v2.0.0)
by Li Lab at The Jackson Laboratory
https://github.com/LabShengLi/nanome
=================================
dsname              : CIEcoli
input               : https://github.com/LabShengLi/nanome/raw/master/test_data/ecoli_ci_test_fast5.tar.gz
genome              : https://storage.googleapis.com/jax-nanopore-01-project-data/nanome-input/ecoli.tar.gz
=================================
[79/843a67] process > EnvCheck (CIEcoli)             [100%] 1 of 1 ✔
[c0/f26448] process > Untar (ecoli_ci_test_fast5.... [100%] 1 of 1 ✔
[29/b92a57] process > Basecall (ecoli_ci_test_fas... [100%] 1 of 1 ✔
[67/ed8b5e] process > QCExport (CIEcoli)             [100%] 1 of 1 ✔
[39/f8cecd] process > Resquiggle (ecoli_ci_test_f... [100%] 1 of 1 ✔
[2e/bd168a] process > Nanopolish (ecoli_ci_test_f... [100%] 1 of 1 ✔
[c9/b4f6c4] process > NplshComb (CIEcoli)            [100%] 1 of 1 ✔
[01/1d622d] process > Megalodon (ecoli_ci_test_fa... [100%] 1 of 1 ✔
[7d/1e4098] process > MgldnComb (CIEcoli)            [100%] 1 of 1 ✔
[b0/c64b3c] process > DeepSignal (ecoli_ci_test_f... [100%] 1 of 1 ✔
[a1/c5c05b] process > DpSigComb (CIEcoli)            [100%] 1 of 1 ✔
[5d/545b2a] process > Report (CIEcoli)               [100%] 1 of 1 ✔
Completed at: 13-May-2022 20:38:30
Duration    : 2m 35s
CPU hours   : 0.2
Succeeded   : 12
```

Results of NANOME pipeline are located at `results/` folder as default:

```
ls results/
CIEcoli-basecallings          CIEcoli_nanome_report.html  README_CIEcoli.txt
CIEcoli-methylation-callings  MultiQC

ls results/CIEcoli-methylation-callings/
Raw_Results-CIEcoli  Read_Level-CIEcoli  Site_Level-CIEcoli  tools_version_table.tsv
```

#### 3.2 NANOME output format

#### Read level output sample is below:
```
zcat results/CIEcoli-methylation-callings/Read_Level-CIEcoli/CIEcoli_NANOME-perRead-score.tsv.gz | head -n 5
ID	Chr	Pos	Strand	Score
21a26cb9-0be2-4670-8694-a3cee91d49b8    NC_000913.3     3503574 -       0.88772296488891
21a26cb9-0be2-4670-8694-a3cee91d49b8    NC_000913.3     3503665 -       2.434816964214898
21a26cb9-0be2-4670-8694-a3cee91d49b8    NC_000913.3     3503670 -       2.434816964214898
21a26cb9-0be2-4670-8694-a3cee91d49b8    NC_000913.3     3503680 -       2.434816964214898
21a26cb9-0be2-4670-8694-a3cee91d49b8    NC_000913.3     3503685 -       2.434816964214898
```
The columns in read level output are:
1. read-id,
1. chromosome,
1. position (1-based),
1. strand,
1. score of log-ratio for probability of 5mC vs. 5C.



#### Site level output sample is below:
```
zcat results/CIEcoli-methylation-callings/Site_Level-CIEcoli/CIEcoli_NANOME-perSite-cov1.sort.bed.gz  | head -n 5
NC_000913.3     3498244 3498245 .       .       +       0.0     1
NC_000913.3     3500081 3500082 .       .       +       0.6666666666666666      3
NC_000913.3     3500082 3500083 .       .       -       1.0     2
NC_000913.3     3500091 3500092 .       .       +       0.6666666666666666      3
NC_000913.3     3500092 3500093 .       .       -       1.0     2
```
The columns in site level output are:
1. chromosome,
1. start (0-based),
1. end (1-based),
1. NA,
1. NA,
1. strand,
1. methylation frequency,
1. coverage.



### 4. More features for NANOME

* Comparisons of state-of-the-art nanopore pipelines
![Figure_pipe_comp](https://raw.githubusercontent.com/LabShengLi/nanome/master/docs/resources/pipeline_comparison.jpg)


* Haplotype-aware consensus methylations. Please check [phasing usage](https://github.com/LabShengLi/nanome/blob/tutorial1/docs/Phasing.md).
![PhasingDemo](https://raw.githubusercontent.com/LabShengLi/nanome/master/docs/resources/nanome3t_5mc_phasing2.png)



### 5. Tutorial video for running NANOME pipeline
[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/TfotM55KTVE/0.jpg)](https://www.youtube.com/watch?v=TfotM55KTVE)

If you have any issues running NANOME, feel free to post at [GitHub](https://github.com/LabShengLi/nanome/issues).


### Reference
1. [https://community.nanoporetech.com](https://community.nanoporetech.com)
2. [https://github.com/nanoporetech/megalodon](https://github.com/nanoporetech/megalodon)
3. [https://github.com/LabShengLi/nanome](https://github.com/LabShengLi/nanome)
4. DNA methylation-calling tools for Oxford Nanopore sequencing: a survey and human epigenome-wide evaluation. Genome Biology 22, 295 (2021). [https://doi.org/10.1186/s13059-021-02510-z](https://doi.org/10.1186/s13059-021-02510-z)
