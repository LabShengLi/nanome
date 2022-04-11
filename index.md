## Tutorial of DNA methylation calling for ONT data
In this totorial, you will learn how to do methylation calling on Oxford Nanopore sequencing data by latest tools. Please create a freshing new folder to execute following commands.

### 1. Software installation
#### 1.1 Install Guppy
Instal the latest version of Guppy:

```
wget https://americas.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_6.1.1_linux64.tar.gz
tar -xzf ont-guppy-cpu_6.1.1_linux64.tar.gz && \
    rm -f ont-guppy-cpu_6.1.1_linux64.tar.gz

ont-guppy-cpu/bin/guppy_basecaller  -v
```

#### 1.2 Install Conda
If you do not have conda, please follow this link [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html) to install conda.

#### 1.3 Install Megalodon
[Megalodon](https://github.com/nanoporetech/megalodon) is a popular and latest ONT developed methylation-calling tool. It can be installed in conda enviroment.

```
conda create --name megalodon python=3.9
conda activate megalodon

pip install megalodon
megalodon -v
```

### 2. Data preparation
In this section, you will prepare the FAST5 files and reference genome for input data.

#### 2.1 Download FAST5 files

```
wget https://github.com/TheJacksonLaboratory/nanome/raw/master/test_data/ecoli_ci_test_fast5.tar.gz
tar -xzf ecoli_ci_test_fast5.tar.gz && \
    rm -f ecoli_ci_test_fast5.tar.gz
    
ls ecoli_ci_test_fast5/
kelvin_20160617_FN_MN17519_sequencing_run_sample_id_74930_ch138_read698_strand.fast5
kelvin_20160617_FN_MN17519_sequencing_run_sample_id_74930_ch139_read4507_strand.fast5
```

#### 2.2 Download reference genome

```
wget https://storage.googleapis.com/jax-nanopore-01-project-data/nanome-input/ecoli.tar.gz
tar -xzf ecoli.tar.gz && rm ecoli.tar.gz

ls ecoli/
Ecoli_k12_mg1655.fasta      Ecoli_k12_mg1655.fasta.bwt           Ecoli_k12_mg1655.fasta.pac
Ecoli_k12_mg1655.fasta.amb  Ecoli_k12_mg1655.fasta.fai           Ecoli_k12_mg1655.fasta.sa
Ecoli_k12_mg1655.fasta.ann  Ecoli_k12_mg1655.fasta.genome.sizes
```

### 3. 5mC & 5hmC detection by Megalodon
Below is an example command to output basecalls, mappings, and CpG 5mC and 5hmC methylation in both per-read (``mod_mappings``) and aggregated (``mods``) formats on prepared data.

```
megalodon \
    ecoli_ci_test_fast5/ \
    --guppy-config  dna_r9.4.1_450bps_fast.cfg\
    --remora-modified-bases dna_r9.4.1_e8 fast 0.0.0 5hmc_5mc CG 0 \
    --outputs basecalls mappings mod_mappings mods per_read_mods \
    --guppy-server-path  ./ont-guppy-cpu/bin/guppy_basecall_server\
    --reference ./ecoli/Ecoli_k12_mg1655.fasta \
    --output-directory methcall_ecoli_data\
    --write-mods-text --overwrite 

ls methcall_ecoli_data/
basecalls.fastq  mappings.summary.txt     per_read_modified_base_calls.db
guppy_log        modified_bases.5hmC.bed  per_read_modified_base_calls.txt
log.txt          modified_bases.5mC.bed   sequencing_summary.txt
mappings.bam     mod_mappings.bam
```

### 4. Consensus 5mC detection by NANOME Pipeline
We developed NANOME, the first Nextflow based container environment (Docker and Singularity) for consensus DNA methylation detection using XGBoost, a gradient boosting algorithm for nanopore long-read sequencing. The consensus outputs can obtain more accurate performance and more comprehensive CpG coverage.

Install Nextflow:

```
conda activate megalodon
# Install nextflow
conda install -c conda-forge -c bioconda nextflow
nextflow -v
```

> :blush: **If you are HPC users, enter into an interactive node before running pipeline.**

Enter an interactive node with 8 cpus for parallelly job running (HPC users only):

```
srun --pty -q batch --time=08:00:00 --mem=25G -n 8  bash
```

Run Nanome consensus pipeline for 5mC detection:

```
module load singularity
mkdir nanome
cd nanome
nextflow run TheJacksonLaboratory/nanome\
    -profile singularity\
    --dsname CIEcoli\
    --input  https://github.com/TheJacksonLaboratory/nanome/raw/master/test_data/ecoli_ci_test_fast5.tar.gz\
    --genome https://storage.googleapis.com/jax-nanopore-01-project-data/nanome-input/ecoli.tar.gz
```

`-profile` is the bundle of paramters used for singularity, `--dsname` is the dataset name, `--input` is the FAST5 input files, and `--genome` is the reference genome file.

The output of NANOME pipeline can be followins:

```
[b2/6dbf5d] process > EnvCheck (CIEcoli)                   [100%] 1 of 1 ✔
[12/5b2e1f] process > Untar (ecoli_ci_test_fast5.tar)      [100%] 1 of 1 ✔
[4f/21a254] process > Basecall (ecoli_ci_test_fast5.tar)   [100%] 1 of 1 ✔
[70/a5ac63] process > QCExport (CIEcoli)                   [100%] 1 of 1 ✔
[cb/2eff17] process > Resquiggle (ecoli_ci_test_fast5.tar) [100%] 1 of 1 ✔
[dc/1747d1] process > Nanopolish (ecoli_ci_test_fast5.tar) [100%] 1 of 1 ✔
[74/8add8f] process > NplshComb (CIEcoli)                  [100%] 1 of 1 ✔
[d1/579953] process > Megalodon (ecoli_ci_test_fast5.tar)  [100%] 1 of 1 ✔
[0e/644818] process > MgldnComb (CIEcoli)                  [100%] 1 of 1 ✔
[07/4265fe] process > DeepSignal (ecoli_ci_test_fast5.tar) [100%] 1 of 1 ✔
[db/653ac8] process > DpSigComb (CIEcoli)                  [100%] 1 of 1 ✔
[0f/7a33dd] process > Guppy (ecoli_ci_test_fast5.tar)      [100%] 1 of 1 ✔
[88/172383] process > GuppyComb (CIEcoli)                  [100%] 1 of 1 ✔
[fc/1e66bb] process > Report (CIEcoli)                     [100%] 1 of 1 ✔
Completed at: 07-Apr-2022 11:13:10
Duration    : 5m 52s
CPU hours   : 0.2
Succeeded   : 14
```

Results of NANOME pipeline are located at `results/` folder as default:

```
ls results/
CIEcoli-basecallings          CIEcoli_nanome_report.html  README_CIEcoli.txt
CIEcoli-methylation-callings  MultiQC

ls results/CIEcoli-methylation-callings/
Raw_Results-CIEcoli  Read_Level-CIEcoli  Site_Level-CIEcoli  tools_version_table.tsv
```

### 5. Tutorial video for running NANOME pipeline
[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/TfotM55KTVE/0.jpg)](https://www.youtube.com/watch?v=TfotM55KTVE)


### Reference
1. https://github.com/nanoporetech/megalodon
2. https://github.com/TheJacksonLaboratory/nanome
3. DNA methylation-calling tools for Oxford Nanopore sequencing: a survey and human epigenome-wide evaluation. Genome Biology 22, 295 (2021). https://doi.org/10.1186/s13059-021-02510-z
