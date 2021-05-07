**This is an explanation of how to use 'nanome' pipeline on raw Fast5 input.**

The inputs of 'nanome' pipeline is a folder/tar/tar.gz or txt file list containing raw signal Fast5 files and a reference genome. We recommend allocate GPU resources to softwares such as Guppy, DeepSignal, DeepMod and Megalodon, in order to optimal running times.

# 1. Running 'nanome' on Fast5 files

The command for running 'nanome' pipeline is to run `./nextflow run https://github.com/liuyangzzu/nano-compare`. `--input` is a compressed file contains Fast5 input file locations, our pipeline support three kinds of inputs: (1) folder, (2) tar/tar.gz file, (3) a txt file `.filelist.txt` contains list of compressed Fast5 files/folders. `--dsname` is output dataset name, `-profile` is the name of execution platform configuration, an example of how to use it is given below.

```angular2html
./nextflow run https://github.com/liuyangzzu/nano-compare \
   --input 'https://github.com/liuyangzzu/nano-compare/raw/6684b71587c19d191ccf87832f25536be10e372f/test_data/demo.fast5.reads.tar.gz' \
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

# 2. Benchmarking experiments
We constructed a list of benchmarking datasets contain Fast5 reads from 800 to 8,000  for NA19240. The datasets can be used upon request. Following command is running 'nanome' pipeline on our benchmarking datasets, reports runing time and memory usage by [Nextflow](https://www.nextflow.io/) utilities. Please refer to the [reports](https://github.com/liuyangzzu/nanome/blob/reproduce-prepare/docs/reports.pdf) and [timeline](https://github.com/liuyangzzu/nanome/blob/reproduce-prepare/docs/timeline.pdf) of benchmarking results on our HPC.

```angular2html
git clone https://github.com/liuyangzzu/nanome.git
cd nanome
./nextflow run main.nf -profile winter  \
	-with-report -with-timeline -with-trace -with-dag \
	-config conf/benchmarking.config
```


We will update more examples here within a short time.