**This is an explanation of how to use 'nanome' pipeline on raw Fast5 input.**

The inputs of 'nanome' pipeline is a folder/tar/tar.gz or txt file list containing raw signal Fast5 files and a reference genome. We recommend allocate GPU resources to softwares such as Guppy, DeepSignal, DeepMod and Megalodon, in order to optimal running times. For METEORE results, it is depends on other tools's read-level outputs (e.g., Megalodon and DeepSignal), and running METEORE program directly on them, detailed please check [METEORE](https://github.com/comprna/METEORE).

# 1. Running 'nanome' for human Nanopore sequencing data

The command for running 'nanome' pipeline is to run `./nextflow run https://github.com/liuyangzzu/nanome`. `--input` is a compressed file contains Fast5 input file locations, our pipeline support three kinds of inputs: (1) folder, (2) tar/tar.gz file, (3) a txt file `.filelist.txt` contains list of compressed Fast5 files/folders. `--dsname` is output dataset name, `--eval true` indicates running evaluation steps,`-profile` is the name of execution platform configuration, an example of our HPC configuration is the server named as [winter](https://github.com/liuyangzzu/nanome/blob/master/nextflow.config#L109), which will include the HPC config file [conf/hpc.config](https://github.com/liuyangzzu/nanome/blob/master/conf/hpc.config). 

By default, we are using hg38 human reference genome, and you can specify reference genome using parameter `--referenceGenome="reference_genome/hg38/hg38.fasta"`. An example of how to use 'nanome' pipeline is given below.

```angular2html
curl -fsSL get.nextflow.io | bash

./nextflow run https://github.com/liuyangzzu/nanome \
   --input 'https://github.com/liuyangzzu/nanome/raw/master/test_data/demo.fast5.reads.tar.gz' \
   --dsname TestData --eval true -profile winter
```

You can also running the pipeline on CloudOS, using following command options.
```angular2html
nextflow run https://github.com/liuyangzzu/nanome
    --config 'conf/google.config'
    --dsname TestData
    --input 'https://github.com/liuyangzzu/nanome/raw/master/test_data/demo.fast5.reads.tar.gz'
```

Command running results is below.

```angular2html
/projects/li-lab/yang/workspace/nano-compare
N E X T F L O W  ~  version 20.10.0
Launching `main.nf` [cranky_mclean] - revision: 39f38efca0
NANOME - NF PIPELINE (v1.0)
by Li Lab at The Jackon Laboratory
http://nanome.jax.org
=================================
dsname          :TestData
input           :inputs/test.demo.filelist.txt
reference_genome    :reference_genome/hg38/hg38.fasta
chromSizesFile      :reference_genome/hg38/hg38.chrom.sizes
runBasecall     :true
runMethcall     :true
evaluation      :false
=================================
[d2/e1e047] process > EnvCheck (EnvCheck)            [100%] 1 of 1 ✔
[67/e8f2d6] process > Untar (demo2.fast5.reads.tar)  [100%] 2 of 2 ✔
[bc/ee9e6f] process > Basecall (demo2.fast5.reads... [100%] 2 of 2 ✔
[3a/99472f] process > QCExport                       [100%] 1 of 1 ✔
[7b/c8a7aa] process > Guppy (demo.fast5.reads.tar)   [100%] 2 of 2 ✔
[54/c65c92] process > GuppyExtract (demo.fast5.re... [100%] 2 of 2 ✔
[60/665d1e] process > Megalodon (demo.fast5.reads... [100%] 2 of 2 ✔
[75/232d0e] process > Resquiggle (demo2.fast5.rea... [100%] 2 of 2 ✔
[7c/6a4152] process > DeepSignal (demo.fast5.read... [100%] 2 of 2 ✔
[98/c7c926] process > Tombo (demo.fast5.reads.tar)   [100%] 2 of 2 ✔
[cc/484339] process > DeepMod (demo.fast5.reads.tar) [100%] 2 of 2 ✔
[60/a7b7e1] process > Nanopolish (demo2.fast5.rea... [100%] 2 of 2 ✔
[7b/2bc531] process > DpSigComb                      [100%] 1 of 1 ✔
[2d/5ca74b] process > TomboComb                      [100%] 1 of 1 ✔
[af/7ca5e5] process > GuppyComb (1)                  [100%] 1 of 1 ✔
[d2/20cfb9] process > MgldnComb                      [100%] 1 of 1 ✔
[ca/910ed5] process > NplshComb                      [100%] 1 of 1 ✔
[1b/1ec959] process > DpmodComb (1)                  [100%] 1 of 1 ✔
Completed at: 16-Aug-2021 13:39:27
Duration    : 41m 42s
CPU hours   : 1.5
Succeeded   : 18
```


All tools's methlation calling and evaluation results will be output to `outputs` folder by default below.

```angular2html
tree outputs/TestData-methylation-callings

outputs/TestData-methylation-callings
├── TestData.deepmod.C.combine.bed.gz
├── TestData.deepsignal.call_mods.combine.tsv.gz
├── TestData.guppy.fast5mod_site_level.combine.tsv.gz
├── TestData.megalodon.per_read.combine.bed.gz
├── TestData.nanopolish.methylation_calls.combine.tsv.gz
└── TestData.tombo.perReadsStats.combined.bed.gz

tree -L 3  outputs/TestData-nanome-analysis/

outputs/TestData-nanome-analysis/
├── MethCorr-TestData_RRBS
│   ├── Meth_corr_plot_data_bgtruth-TestData_RRBS-bsCov1-minToolCov1-baseFormat1.csv
│   ├── Meth_corr_plot_data_joined-TestData_RRBS-bsCov1-minToolCov1-baseFormat1.csv
│   ├── Meth_corr_plot_data-TestData_RRBS-correlation-matrix.xlsx
│   ├── run-results.log
│   ├── TestData_RRBS-summary-bgtruth-tools-bsCov1-minCov1.csv
│   ├── venn.data.TestData_RRBS.TestData.five.tools.cov1.dat
│   └── venn.data.TestData_RRBS.TestData.top3.cov1.dat
└── MethPerf-TestData_RRBS
    ├── performance-results
    │   ├── curve_data
    │   ├── TestData_RRBS.DeepMod.performance.report.csv
    │   ├── TestData_RRBS.DeepSignal.performance.report.csv
    │   ├── TestData_RRBS.Megalodon.performance.report.csv
    │   └── TestData_RRBS.Nanopolish.performance.report.csv
    ├── run-results.log
    ├── TestData_RRBS.hg38_nonsingletons.concordant.bed
    ├── TestData_RRBS.hg38_nonsingletons.discordant.bed
    ├── TestData_RRBS.summary.bsseq.singleton.nonsingleton.cov1.csv
    └── TestData_RRBS.Tools_BGTruth_cov1_Joined.bed
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