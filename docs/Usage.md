**This is an explanation of how to use 'nanome' pipeline on raw Fast5 input.**

The inputs of 'nanome' pipeline is a folder/tar/tar.gz or txt file list containing raw signal Fast5 files and a reference genome. We recommend allocate GPU resources to softwares such as Guppy, DeepSignal, DeepMod and Megalodon, in order to optimal running times.

# 1. Running 'nanome' for human Nanopore sequencing data

The command for running 'nanome' pipeline is to run `./nextflow run https://github.com/liuyangzzu/nanome`. `--input` is a compressed file contains Fast5 input file locations, our pipeline support three kinds of inputs: (1) folder, (2) tar/tar.gz file, (3) a txt file `.filelist.txt` contains list of compressed Fast5 files/folders. `--dsname` is output dataset name, `-profile` is the name of execution platform configuration. 

By default, we are using hg38 human reference genome, and you can specify reference genome using parameter `--referenceGenome="reference_genome/hg38/hg38.fasta"`. An example of how to use 'nanome' pipeline is given below.

```angular2html
curl -fsSL get.nextflow.io | bash

./nextflow run https://github.com/liuyangzzu/nanome \
   --input 'https://github.com/liuyangzzu/nanome/raw/master/test_data/demo.fast5.reads.tar.gz' \
   --dsname TestData -profile winter
```

All tools's methlation calling and evaluation results will be output to `outputs` folder by default below.

```angular2html
tree outputs/TestData-methylation-callings

outputs/TestData-methylation-callings
├── TestData.DeepModC_clusterCpG.combine.bed.gz
├── TestData.DeepSignal.combine.tsv.gz
├── TestData.Megalodon.combine.bed.gz
├── TestData.Nanopolish.combine.tsv.gz
└── TestData.Tombo.combine.bed.gz

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
# 2. Experiment for E. coli data
The 'nanome' pipeline supports 5mC detection by all tools on both human and Escherichia coli data. Note that `--referenceGenome` need to be set as E. coli reference genome such as 'reference_genome/ecoli/Ecoli_k12_mg1655.fasta'. Below is an example of pipeline runing on E. coli data, please refer to the input parameters for pipeline `-config` params [conf/ecoli_demo.config](https://github.com/liuyangzzu/nanome/blob/master/conf/ecoli_demo.config).

```angular2html
git clone https://github.com/liuyangzzu/nanome.git
cd nanome
curl -fsSL get.nextflow.io | bash

./nextflow run main.nf -profile winter \
    -config conf/ecoli_demo.config
```

The outputs of pipeline on E. coli data is results by all tools methylation calling below.

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
We constructed a list of benchmarking datasets contain Fast5 reads from 800 to 8,000  for NA19240. The datasets can be used upon request. Following command is running 'nanome' pipeline on our benchmarking datasets. 

```angular2html
git clone https://github.com/liuyangzzu/nanome.git
cd nanome

nextflow run main.nf -profile winter  \
	-with-report -with-timeline -with-trace -with-dag \
	-config conf/benchmarking.config
```

Reports runing time and memory usage by [Nextflow](https://www.nextflow.io/) utilities. Please refer to the [reports](https://github.com/liuyangzzu/nanome/blob/master/docs/reports.pdf) and [timeline](https://github.com/liuyangzzu/nanome/blob/master/docs/timeline.pdf) of benchmarking results on our HPC.

We will update more examples here within a short time.