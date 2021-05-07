# DNA methylation calling tools for Oxford Nanopore sequencing: a survey and human epigenome-wide evaluation
## -- The 'nanome' pipeline for evaluation of DNA methylation calling tools for Oxford Nanopore sequencing 

## Methodology of nanome pipeline

**Background:** Nanopore long-read sequencing technology greatly expands the capacity of long-range single-molecule DNA-modification detection. A growing number of analytical tools have been actively developed to detect DNA methylation from Nanopore sequencing reads. Here, we examine the performance of different methylation calling tools to provide a systematic evaluation to guide practitioners for human epigenome-wide research.


![Figure1](https://github.com/liuyangzzu/nanome/blob/reproduce-prepare/docs/Fig1.jpg)
**Fig 1.**  **(A)** Timeline of publication and technological developments of Oxford Nanopore Technologies (ONT) methylation calling tools to detect DNA cytosine modifications. **(B)** Performance evaluation on 5mC/5C prediction of methylation calling tools with Nanopore sequencing.


![FigureS1](https://github.com/liuyangzzu/nanome/blob/reproduce-prepare/docs/fig.s2.workflow.jpg)
**Fig S1. Workflow for 5-methylcytosine (5mC) detection for Nanopore sequencing.** The pipeline has four steps: (1) Nanopore sequencing. (2) Base calling, which requires raw signals and reference genome as input to perform base calling by Albacore or Guppy. (3) Alignment to the genome by direct mapping with miniMap2 and re-squiggle with Tombo (optional). (4) Methylation calling. Here we compare five methylation calling tools: Nanopolish, Megalodon, DeepSignal, Tombo, and DeepMod to detect cytosine status in CpG context.


**Results:** We compare five analytic frameworks for detecting DNA modification from Nanopore long-read sequencing data. We evaluate the association between genomic context, CpG methylation-detection accuracy, CpG sites coverage, and running time using Nanopore sequencing data from natural human DNA. Furthermore, we provide an online DNA methylation database (https://nanome.jax.org) with which to display genomic regions that exhibit differences in DNA-modification detection power among different methylation calling algorithms for nanopore sequencing data.


![Figure4](https://github.com/liuyangzzu/nanome/blob/reproduce-prepare/docs/Fig4.jpg)
**Fig 4. Comparison of Nanopore methylation calling tools for the detection of CpG methylation on four real world data sets in biologically relevant genomic contexts.**


![Figure5](https://github.com/liuyangzzu/nanome/blob/reproduce-prepare/docs/Fig5.jpg)
**Fig 5. Comparison of Pearson correlation of methylation call tools across all CpG sites.** **(A)** Correlation plot showing Pearson correlation of each methylation calling tool with BS- Seq on NA19240. **(B-C)** Relationship between CpG methylation percentage and distance to annotated TSS in (B) NA19240 and (C) APL. **(D)** Relationship between CpG methylation percentage and distance to annotated CTCF binding peaks in NA19240. Distances are binned into (B, C) 50-bp, and (D) 100-bp windows. Negative distances are upstream and positive distances are downstream of the (B-C) TSS and CTCF binding peaks (D). 


![Figure6](https://github.com/liuyangzzu/nanome/blob/reproduce-prepare/docs/Fig6.jpg)
**Fig 6. CpG sites detected by methylation calling tools using NA19240.** UpSet diagram shown at the lower left is for CpG sites detected by all methylation calling tools. Venn diagram shown at the upper right is for CpG sites detected by Top3 performance methylation calling tools (Nanopolish, Megalodon and DeepSignal). For each methylation calling tool, only CpG sites covered >= 3 reads are considered.


![Figure7](https://github.com/liuyangzzu/nanome/blob/reproduce-prepare/docs/Fig7.jpg)
**Fig 7. CPU utilized time and memory usage for each methylation calling tool on each dataset.** All tools were compared on the same computer clusters: 32 cores, 2.6GHz HP Proliant SL Series CPU, 300 GB RAM, NVIDIA Tesla P100 Data Center and 1 TB Data Direct Networks Gridscalar GS7k GPFS storage appliance. The HPC platform software and hardware specifications are: slurm manager version: 19.05.5, CPU: Intel(R) Xeon(R) Gold 6136 CPU @ 3.00GHz, GPU: Tesla V100-SXM2-32GB.

**Benchmarking reports using [Nextflow](https://www.nextflow.io/)**

We construct a benchmarking data contains reads from 800 to 8,000 reads for NA19240, and monitoring job running timeline and resource usage on our HPC, reports generated nextflow are: [Report](https://github.com/liuyangzzu/nanome/blob/reproduce-prepare/docs/reports.pdf)  and [Timeline](https://github.com/liuyangzzu/nanome/blob/reproduce-prepare/docs/timeline.pdf). 

**Conclusions:** Our study is the first benchmark of computational methods for mammalian whole genome DNA-modification detection in Nanopore sequencing. We provide a broad foundation for cross-platform standardization, and an evaluation of analytical tools designed for genome-scale modified-base detection using Nanopore sequencing.


## Revision History

For release history, please visit [here](https://github.com/liuyangzzu/nanome/releases). For details, please go [here](https://github.com/liuyangzzu/nanome/blob/reproduce-prepare/README.md).

## Contact

If you have any questions/issues/bugs, please post them on [GitHub](https://github.com/liuyangzzu/nanome/issues). We will regularly update the Github to support more methylation calling tools.

## Reference
**Detailed results can be found in our preprint. Please cite our preprint below if you are interested in our pipeline:**

Yang Liu, Wojciech Rosikiewicz, Ziwei Pan, Nathaniel Jillette, Aziz Taghbalout, Jonathan Foox, Christopher Mason, Martin Carroll, Albert Cheng, Sheng Li. DNA methylation calling tools for Oxford Nanopore sequencing: a survey and human epigenome-wide evaluation. bioRxiv, 2020. Online at https://doi.org/10.1101/2021.05.05.442849.

