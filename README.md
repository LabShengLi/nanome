# DNA methylation calling tools for Oxford Nanopore sequencing: a survey and human epigenome-wide evaluation
## -- The 'nanome' pipeline for evaluation of DNA methylation calling tools for Oxford Nanopore sequencing 

## Methodology of nanome pipeline

**Background:** Nanopore long-read sequencing technology greatly expands the capacity of long-range single-molecule DNA-modification detection. A growing number of analytical tools have been actively developed to detect DNA methylation from Nanopore sequencing reads. Here, we examine the performance of different methylation calling tools to provide a systematic evaluation to guide practitioners for human epigenome-wide research.


![Figure1](https://github.com/liuyangzzu/nanome/blob/reproduce-prepare/docs/Fig1.jpg)
**Fig 1.**  **(A)** Timeline of publication and technological developments of Oxford Nanopore Technologies (ONT) methylation calling tools to detect DNA cytosine modifications. **(B)** Performance evaluation on 5mC/5C prediction of methylation calling tools with Nanopore sequencing.


![FigureS1](https://github.com/liuyangzzu/nanome/blob/reproduce-prepare/docs/fig.s2.workflow.jpg)
**Fig S1. Workflow for 5-methylcytosine (5mC) detection for Nanopore sequencing.** The pipeline has four steps: (1) Nanopore sequencing. (2) Base calling, which requires raw signals and reference genome as input to perform base calling by Albacore or Guppy. (3) Alignment to the genome by direct mapping with miniMap2 and re-squiggle with Tombo (optional). (4) Methylation calling. Here we compare five methylation calling tools: Nanopolish, Megalodon, DeepSignal, Tombo, and DeepMod to detect cytosine status in CpG context.


**Results:** We compare five analytic frameworks for detecting DNA modification from Nanopore long-read sequencing data. We evaluate the association between genomic context, CpG methylation-detection accuracy, CpG sites coverage, and running time using Nanopore sequencing data from natural human DNA. Furthermore, we provide an online DNA methylation database (https://nanome.jax.org) with which to display genomic regions that exhibit differences in DNA-modification detection power among different methylation calling algorithms for nanopore sequencing data.

**Conclusions:** Our study is the first benchmark of computational methods for mammalian whole genome DNA-modification detection in Nanopore sequencing. We provide a broad foundation for cross-platform standardization, and an evaluation of analytical tools designed for genome-scale modified-base detection using Nanopore sequencing.


## Revision History

For release history, please visit [here](https://github.com/liuyangzzu/nanome/releases). For details, please go [here](https://github.com/liuyangzzu/nanome/blob/reproduce-prepare/README.md).

## Contact

If you have any questions/issues/bugs, please post them on [GitHub](https://github.com/liuyangzzu/nanome/issues). 

## Reference
**Please cite the publication below if you use our pipelines:**

Yang Liu, Wojciech Rosikiewicz, Ziwei Pan, Nathaniel Jillette, Aziz Taghbalout, Jonathan Foox, Christopher Mason, Martin Carroll, Albert Cheng, Sheng Li. DNA methylation calling tools for Oxford Nanopore sequencing: a survey and human epigenome-wide evaluation. bioRxiv, 2020. Online at https://doi.org/10.1101/2021.05.05.442849.

