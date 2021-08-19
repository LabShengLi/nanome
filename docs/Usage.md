**This is an explanation of how to use 'nanome' pipeline on raw Fast5 input.**

The inputs of 'nanome' pipeline is a folder/tar/tar.gz or txt file list containing raw signal Fast5 files and a reference genome. We recommend allocate GPU resources to softwares such as Guppy, DeepSignal, DeepMod and Megalodon, in order to optimal running times. For getting METEORE result, it is depends on other tools's read-level outputs (e.g., Megalodon and DeepSignal), and running METEORE program directly on them, detailed please check [METEORE](https://github.com/comprna/METEORE).

# 1. Running 'nanome' for human Nanopore sequencing data

The command for running 'nanome' pipeline is to run `./nextflow run https://github.com/liuyangzzu/nanome`. `--input` is a compressed file contains Fast5 input file locations, our pipeline support three kinds of inputs: (1) folder, (2) tar/tar.gz file, (3) a txt file `.filelist.txt` contains list of compressed Fast5 files/folders. `--dsname` is output dataset name, `--eval true` indicates running evaluation steps,`-profile` is the name of execution platform configuration, an example of our HPC configuration is the server named as [winter](https://github.com/liuyangzzu/nanome/blob/master/nextflow.config#L109), which will include the HPC config file [conf/hpc.config](https://github.com/liuyangzzu/nanome/blob/master/conf/hpc.config). 

By default, we are using hg38 human reference genome, and you can specify reference genome using parameter `--referenceGenome="reference_genome/hg38/hg38.fasta"`. An example of how to use 'nanome' pipeline is given below.

```angular2html
# Get nextflow executable file
curl -fsSL get.nextflow.io | bash

# Get nanome singularity
singularity pull nanome_v1.4.sif docker://quay.io/liuyangzzu/nanome:v1.4

# Run nanome pipeline on project directory
./nextflow run main.nf \
    -profile winter2 \
    -with-report -with-timeline -with-trace -with-dag \
    -with-singularity nanome_v1.4.sif \
    --dsname TestData
    --input https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt

# Running nanome pipeline directly from github
nextflow run https://github.com/liuyangzzu/nanome.git \
    -profile winter2 \
    -with-singularity nanome_v1.4.sif \
    --dsname TestData \
    --input https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt
```

You can also running the pipeline on CloudOS, using following command options.
```angular2html
nextflow run https://github.com/liuyangzzu/nanome.git \
    --config 'conf/google.config' \
    --dsname TestData \
    --input https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt
```

Pipeline running results is below, please also check the pipeline output directory tree for [outputs](https://github.com/liuyangzzu/nanome/blob/doc-task/docs/outputs.tree.txt) and [work](https://github.com/liuyangzzu/nanome/blob/doc-task/docs/work.tree.txt). It can also generates [timeline](https://github.com/liuyangzzu/nanome/blob/doc-task/docs/resources/timeline_demo.pdf), [report](https://github.com/liuyangzzu/nanome/blob/doc-task/docs/resources/report_demo.pdf) and [resource usage](https://github.com/liuyangzzu/nanome/blob/doc-task/docs/resources/trace_demo.txt.tsv).

```angular2html
N E X T F L O W  ~  version 20.10.0
Launching `main.nf` [hungry_bartik] - revision: 12c6a8e37a
NANOME - NF PIPELINE (v1.0)
by Li Lab at The Jackson Laboratory
https://nanome.jax.org
=================================
dsname          :TestData
input           :https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt
reference_genome    :reference_genome/hg38/hg38.fasta
chromSizesFile      :reference_genome/hg38/hg38.chrom.sizes
runBasecall     :true
runMethcall     :true
evaluation      :true
=================================
executor >  slurm (30)
[a8/145382] process > EnvCheck (EnvCheck)            [100%] 1 of 1 ✔
[94/d0e364] process > Untar (demo.fast5.reads.tar)   [100%] 2 of 2 ✔
[e4/a6c42f] process > Basecall (demo.fast5.reads.... [100%] 2 of 2 ✔
[17/87d46b] process > QCExport                       [100%] 1 of 1 ✔
[f0/8b1891] process > Guppy (demo.fast5.reads.tar)   [100%] 2 of 2 ✔
[d4/2f7794] process > Megalodon (demo.fast5.reads... [100%] 2 of 2 ✔
[3f/336996] process > Resquiggle (demo.fast5.read... [100%] 2 of 2 ✔
[78/76a38d] process > DeepSignal (demo.fast5.read... [100%] 2 of 2 ✔
[b4/cd7309] process > Tombo (demo.fast5.reads.tar)   [100%] 2 of 2 ✔
[21/92aae7] process > DeepMod (demo.fast5.reads.tar) [100%] 2 of 2 ✔
[bd/86ab7c] process > Nanopolish (demo.fast5.read... [100%] 2 of 2 ✔
[aa/de3cba] process > DpSigComb                      [100%] 1 of 1 ✔
[67/f0ba5b] process > TomboComb                      [100%] 1 of 1 ✔
[06/911c00] process > GuppyComb (1)                  [100%] 1 of 1 ✔
[71/7a6fe9] process > MgldnComb                      [100%] 1 of 1 ✔
[4c/70dbe9] process > NplshComb                      [100%] 1 of 1 ✔
[ae/17edb3] process > DpmodComb (1)                  [100%] 1 of 1 ✔
[cf/706926] process > METEORE (1)                    [100%] 1 of 1 ✔
[a5/732325] process > SiteLevelUnify (1)             [100%] 1 of 1 ✔
[77/f19512] process > ReadLevelPerf (1)              [100%] 1 of 1 ✔
[37/3cdaa5] process > SiteLevelCorr (1)              [100%] 1 of 1 ✔
Completed at: 18-Aug-2021 16:04:14
Duration    : 1h 46s
CPU hours   : 1.9
Succeeded   : 30
```


All tools's methlation calling and evaluation results will be output to `outputs` folder by default below.

```angular2html
tree outputs/TestData-methylation-callings/

outputs/TestData-methylation-callings/
├── TestData.deepmod.C_clusterCpG_per_site.combine.bed.gz
├── TestData.deepsignal.per_read.combine.tsv.gz
├── TestData.guppy.fast5mod_per_site.combine.tsv.gz
├── TestData.megalodon.per_read.combine.bed.gz
├── TestData.meteore.megalodon_deepsignal_optimized_model_per_read.combine.tsv.gz
├── TestData.nanopolish.per_read.combine.tsv.gz
└── TestData.tombo.per_read.combine.bed.gz

tree outputs/TestData-qc-report

outputs/TestData-qc-report
├── TestData.coverage.negativestrand.bed.gz
├── TestData.coverage.positivestrand.bed.gz
└── TestData-qc-report.tar.gz

tree outputs/nanome-analysis-TestData/ -L 2

outputs/nanome-analysis-TestData/
├── MethCorr-TestData_RRBS
│   ├── DONE.txt
│   ├── Meth_corr_plot_data_bgtruth-TestData_RRBS-bsCov1-minToolCov1-baseFormat1.csv.gz
│   ├── Meth_corr_plot_data_joined-TestData_RRBS-bsCov1-minToolCov1-baseFormat1.csv.gz
│   ├── Meth_corr_plot_data-TestData_RRBS-correlation-matrix-toolcov1-bsseqcov1.xlsx
│   ├── run-results.log
│   ├── TestData.tools.cov1.join.with.bsseq.cov1.site.level.report.csv
│   └── venn_data
├── MethPerf-TestData_RRBS
│   ├── DONE.txt
│   ├── performance-results
│   ├── run-results.log
│   ├── TestData_RRBS.hg38_nonsingletons.concordant.bed.gz
│   ├── TestData_RRBS.hg38_nonsingletons.discordant.bed.gz
│   ├── TestData_RRBS.summary.bsseq.cov1.joined.tools.singleton.nonsingleton.table.like.s2.csv
│   ├── TestData_RRBS.summary.bsseq.singleton.nonsingleton.cov1.csv
│   ├── TestData_RRBS.Tools_BGTruth_cov1_Joined_baseFormat1.bed.gz
│   └── TestData.tools.cov1.join.with.bsseq.cov1.read.level.report.xlsx
├── Read_Level-TestData
│   ├── TestData_DeepSignal-METEORE-perRead-score.tsv.gz
│   ├── TestData_Guppy-METEORE-perRead-score.tsv.gz
│   ├── TestData_Megalodon-METEORE-perRead-score.tsv.gz
│   ├── TestData_Nanopolish-METEORE-perRead-score.tsv.gz
│   └── TestData_Tombo-METEORE-perRead-score.tsv.gz
└── Site_Level-TestData
    ├── Site_Level-TestData.tss.DeepMod.cov1.bed.gz
    ├── Site_Level-TestData.tss.DeepSignal.cov1.bed.gz
    ├── Site_Level-TestData.tss.Guppy.cov1.bed.gz
    ├── Site_Level-TestData.tss.Megalodon.cov1.bed.gz
    ├── Site_Level-TestData.tss.METEORE.cov1.bed.gz
    ├── Site_Level-TestData.tss.Nanopolish.cov1.bed.gz
    └── Site_Level-TestData.tss.Tombo.cov1.bed.gz
```

We also support input as a file list if input is suffix like `.filelist.txt`, an example input is [inputs/test.demo.filelist.txt](https://github.com/liuyangzzu/nanome/blob/master/inputs/test.demo.filelist.txt).
```angular2html
./nextflow run main.nf -profile winter --input inputs/test.demo.filelist.txt
```

# 2. Experiment for E. coli data
The 'nanome' pipeline supports 5mC detection by all tools on both human and Escherichia coli data. Note that `--referenceGenome` need to be set as E. coli reference genome such as 'reference_genome/ecoli/Ecoli_k12_mg1655.fasta'. Below is an example of pipeline runing on E. coli data, please refer to the input parameters for pipeline `-config` params [conf/ecoli_demo.config](https://github.com/liuyangzzu/nanome/blob/master/conf/ecoli_demo.config).

```angular2html
git clone https://github.com/liuyangzzu/nanome.git
cd nanome
curl -fsSL get.nextflow.io | bash

./nextflow run main.nf -profile winter \
    -config conf/ecoli_demo.config
```
Command running results is below.

```angular2html
N E X T F L O W  ~  version 20.10.0
Launching `main.nf` [friendly_leakey] - revision: eafc216253
NANOME - NF PIPELINE (v1.0)
by Li Lab at The Jackon Laboratory
http://nanome.jax.org
=================================
dsname              :EcoliDemo
input               :/projects/li-lab/Nanopore_compare/suppdata/ecoli-sanity-check/ecoli_meteore.tar.gz
reference_genome    :reference_genome/ecoli/Ecoli_k12_mg1655.fasta
chromSizesFile      :reference_genome/ecoli/Ecoli_k12_mg1655.fasta.genome.sizes
runBasecall         :true
runMethcall         :true
eval                :false
=================================
executor >  slurm (14)
[3f/f6f1d7] process > EnvCheck (EnvCheck)                                  [100%] 1 of 1 ✔
[f8/ce877c] process > Basecall (ecoli_meteore)                             [100%] 1 of 1 ✔
[d5/f84584] process > QCExport                                             [100%] 1 of 1 ✔
[21/74324c] process > Resquiggle (ecoli_meteore_basecalled)                [100%] 1 of 1 ✔
[2f/f36517] process > DeepSignal (ecoli_meteore_basecalled_resquiggle_dir) [100%] 1 of 1 ✔
[ba/166a6f] process > Tombo (ecoli_meteore_basecalled_resquiggle_dir)      [100%] 1 of 1 ✔
[5b/c3613f] process > Megalodon (ecoli_meteore)                            [100%] 1 of 1 ✔
[2a/f80bc3] process > DeepMod (ecoli_meteore_basecalled)                   [100%] 1 of 1 ✔
[29/d9b126] process > Nanopolish (ecoli_meteore_basecalled)                [100%] 1 of 1 ✔
[2b/18376c] process > DpSigComb                                            [100%] 1 of 1 ✔
[20/ca5232] process > TomboComb                                            [100%] 1 of 1 ✔
[c7/4c2e8d] process > MgldnComb                                            [100%] 1 of 1 ✔
[4c/50e8a1] process > NplshComb                                            [100%] 1 of 1 ✔
[c5/ad8ca7] process > DpmodComb                                            [100%] 1 of 1 ✔
Completed at: 08-May-2021 16:53:56
Duration    : 7m 35s
CPU hours   : 0.4
Succeeded   : 14
```

The output files of pipeline on E. coli data by all tools are below.

```angular2html
tree outputs/EcoliDemo-methylation-callings

outputs/EcoliDemo-methylation-callings
├── EcoliDemo.DeepModC.combine.bed.gz
├── EcoliDemo.DeepSignal.combine.tsv.gz
├── EcoliDemo.Megalodon.combine.bed.gz
├── EcoliDemo.Nanopolish.combine.tsv.gz
└── EcoliDemo.Tombo.combine.bed.gz
```


# 3. Benchmarking experiment
We constructed a list of benchmarking datasets that contain Fast5 reads from 800 to 8,000  for NA19240. The datasets can be got by users upon request. Following command is running 'nanome' pipeline on our benchmarking datasets, please refer to the input parameters for pipeline `-config` params [conf/benchmarking.config](https://github.com/liuyangzzu/nanome/blob/master/conf/benchmarking.config).

```angular2html
git clone https://github.com/liuyangzzu/nanome.git
cd nanome

nextflow run main.nf -profile winter  \
	-with-report -with-timeline -with-trace -with-dag \
	-config conf/benchmarking.config
```

Resource usage are reported by [Nextflow](https://www.nextflow.io/) workflow reporting utilities. Please refer to the [Trace file](https://github.com/liuyangzzu/nanome/blob/master/docs/nanome.pipeline_trace.tsv), [Report](https://github.com/liuyangzzu/nanome/blob/master/docs/reports2.pdf) and [Timeline](https://github.com/liuyangzzu/nanome/blob/master/docs/timeline.pdf) of benchmarking results on our HPC.

# 4. Running pipeline on cloud computing platform

Our Nextflow pipeline can running on CloudOS. The CloudOS will use a default Docker image. Below is an example.

```angular2html
git clone https://github.com/liuyangzzu/nanome.git
cd nanome
curl -s https://get.nextflow.io | bash

./nextflow run main.nf -profile gls -w gs://jax-nanopore-01-project-data/nanome-work-test --outputDir gs://jax-nanopore-01-project-data/nanome-outputs
```

The `jax-nanopore-01-project-data` is a sample of **Data Bucket** name that you can access on google cloud. `-w` is pipeline output working directory, `--outputDir` is methylation calling and evaluation results directory.

For more detail of using cloud computing, please check [Cloud computing usage](https://github.com/liuyangzzu/nanome/blob/master/docs/CloudComputing.md).

We will update more examples here within a short time.