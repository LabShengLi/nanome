**This is an explanation of how to use script to perform read-level and site-level performance evaluation.**

User needs to provide ONT tools's methylation-calling raw outputs and BS-seq data before he can perform evaluations. The genome-annotation files may needed if the performance at specific genomic regions is interested.

You can install NANOME evaluation script from [PyPI](https://pypi.org/project/nanome-jax/):
```angular2html
pip install nanome-jax
```

# 1. Read-level performance evaluation
## Script for read-level evaluation
The script `read_level_eval.py` is desinged for general purpose of read-level performance evaluation for all kinds of tools.

```angular2html
read_level_eval.py -v

read_level_eval (NANOME) v1.3.2
```

## Sample usage
User needs to specify methylation-calling raw results of each tool and BS-seq replicates to execute the read-level evaluation script.
```angular2html
read_level_eval.py \
    --calls\
        Nanopolish:[Nanopolish-calls]\
        Megalodon:[Megalodon-calls]\
        DeepSignal:[DeepSignal-calls]\
    --bgtruth bismark:[BS-seq-rep1];[BS-seq-rep2]\
    --runid MethPerf-TestData_RRBS_2Reps\
    --dsname TestData\
    --min-bgtruth-cov 1\
    --report-joined\
    --genome-annotation [genome-annotation]\
    -o .
```
Sample results can be found at [read-level outputs](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/read_level_output.txt).

## Command options
```angular2html
read_level_eval.py -h

usage: read_level_eval (NANOME) [-h] [-v] [--dsname DSNAME] --runid RUNID
                                --calls CALLS [CALLS ...] [--bgtruth BGTRUTH]
                                [--min-bgtruth-cov MIN_BGTRUTH_COV]
                                [--processors PROCESSORS] [--report-joined]
                                [--chrSet CHRSET [CHRSET ...]] [-o O]
                                [--enable-cache] [--using-cache]
                                [--distribution] [--mpi] [--analysis ANALYSIS]
                                [--bedtools-tmp BEDTOOLS_TMP]
                                [--genome-annotation GENOME_ANNOTATION]
                                [--reference-genome REFERENCE_GENOME]
                                [--save-curve-data] [--large-mem]
                                [--disable-bed-check] [--debug]

Read-level performance evaluation in nanome paper

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --dsname DSNAME       dataset name
  --runid RUNID         running prefix
  --calls CALLS [CALLS ...]
                        all ONT call results <tool-name>:<file-name> seperated
                        by space
  --bgtruth BGTRUTH     background truth file <encode-type>:<file-name>;<file-
                        name>
  --min-bgtruth-cov MIN_BGTRUTH_COV
                        min bg-truth coverage cutoff
  --processors PROCESSORS
                        multi-processors
  --report-joined       true if report on only joined sets
  --chrSet CHRSET [CHRSET ...]
                        chromosome list
  -o O                  output dir
  --enable-cache        if enable cache functions
  --using-cache         if use cache files
  --distribution        report singleton/nonsingleton distributions
  --mpi                 if using multi-processing for evaluation
  --analysis ANALYSIS   special analysis specifications for ecoli
  --bedtools-tmp BEDTOOLS_TMP
                        bedtools temp dir
  --genome-annotation GENOME_ANNOTATION
                        genome annotation dir
  --reference-genome REFERENCE_GENOME
                        genome annotation dir
  --save-curve-data     if save pred/truth points for curve plot
  --large-mem           if using large memory (>100GB)
  --disable-bed-check   if disable checking the 0/1 base format for genome
                        annotations
  --debug               if output debug info
```

# 2. Site-level performance evaluation

## Script for site-level evaluation
The script `site_level_eval.py` is desinged for general purpose of site-level performance evaluation for all kinds of tools.

```angular2html
site_level_eval.py -v

site_level_eval (NANOME) v1.3.2
```


## Sample usage
User needs to specify methylation-calling raw results of each tool and BS-seq replicates to execute the site-level evaluation script.

```angular2html
site_level_eval.py \
    --calls\
        Nanopolish:[Nanopolish-calls]\
        Megalodon:[Megalodon-calls]\
        DeepSignal:[DeepSignal-calls]\
    --bgtruth bismark:[BS-seq-rep1];[BS-seq-rep2]\
    --runid MethCorr-TestData_RRBS_2Reps\
    --dsname TestData\
    --min-bgtruth-cov 3 --toolcov-cutoff 1\
    --genome-annotation [genome-annotation]\
    --beddir .\
    --gen-venn --summary-coverage\
    -o . 
```

Sample results can be found at [site-level outputs](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/resources/site_level_output.txt).

## Command options
```angular2html
site_level_eval.py -h

usage: site_level_eval (NANOME) [-h] [-v] --dsname DSNAME --runid RUNID
                                --calls CALLS [CALLS ...] --bgtruth BGTRUTH
                                [--beddir BEDDIR]
                                [--min-bgtruth-cov MIN_BGTRUTH_COV]
                                [--toolcov-cutoff TOOLCOV_CUTOFF] [--sep SEP]
                                [--processors PROCESSORS] [-o O] [--gen-venn]
                                [--summary-coverage] [--region-coe-report]
                                [--enable-cache] [--using-cache] [--plot]
                                [--bedtools-tmp BEDTOOLS_TMP]
                                [--genome-annotation GENOME_ANNOTATION]
                                [--large-mem] [--disable-bed-check] [--debug]

Site-level correlation analysis in nanome paper

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --dsname DSNAME       dataset name
  --runid RUNID         running prefix
  --calls CALLS [CALLS ...]
                        all ONT call results <tool-name>:<file-name> seperated
                        by spaces
  --bgtruth BGTRUTH     background truth file <encode-type>:<file-
                        name1>;<file-name1>
  --beddir BEDDIR       base dir for concordant/discordant bed files
  --min-bgtruth-cov MIN_BGTRUTH_COV
                        cutoff for coverage in bg-truth
  --toolcov-cutoff TOOLCOV_CUTOFF
                        cutoff for coverage in nanopore tools
  --sep SEP             seperator for output csv file
  --processors PROCESSORS
                        running processors
  -o O                  output dir
  --gen-venn            generate CpGs for venn data analysis
  --summary-coverage    generate summary table for coverage at each region
  --region-coe-report   report table for PCC value at each region
  --enable-cache        if enable cache functions
  --using-cache         if use cache files
  --plot                plot the correlation matrix figure
  --bedtools-tmp BEDTOOLS_TMP
                        bedtools temp dir
  --genome-annotation GENOME_ANNOTATION
                        genome annotation dir
  --large-mem           if using large memory (>100GB)
  --disable-bed-check   if disable checking the 0/1 base format for genome
                        annotations
  --debug               if output debug info
```

# 3. Unified read/site level format conversion


## Script for unified format conversion
The script `tss_eval.py` is desinged for general purpose of converting raw results of all kinds of tools into a unified format.

```angular2html
tss_eval.py -v

tss_eval (NANOME) v1.3.2
```

## Sample usage for read-level format unification
The option `--output-unified-format` is used for generate unified read-level format for tools.
```angular2html
tss_eval.py \
    --calls \
        [tool-name]:[tool-call-filename] \
    --runid Read_Level-[dataset-name] \
    --dsname [dataset-name]\
    --output-unified-format\
    -o [output-dir]
```

The unified format of read-level output is below, the `Pos` is 1-base position:
```angular2html
ID	Chr	Pos	Strand	Score
2414d963-488b-4987-9b28-6c5f2af76e4e	chr1	125179705	+	2.3647238093176295
2414d963-488b-4987-9b28-6c5f2af76e4e	chr1	125179715	+	3.3257023748112737
2414d963-488b-4987-9b28-6c5f2af76e4e	chr1	125179724	+	2.5346507818507837
2414d963-488b-4987-9b28-6c5f2af76e4e	chr1	125179734	+	0.42996208256507207
2414d963-488b-4987-9b28-6c5f2af76e4e	chr1	125179754	+	0.8669499128929393
```

## Sample usage for site-level format unification
Below will generate the site level unified format for tools.
```angular2html
tss_eval.py \
    --calls \
        [tool-name]:[tool-call-filename]\
    --runid Site_Level-[dataset-name]\
    --dsname [dataset-name]\
    -o [output-dir]

```

The unified format of site-level output is below, the first three columns are chr, 0-base start, 1-base end, and the last three columns are strand, methylation frequency and coverage.
```angular2html
chr1	11343977	11343978	.	.	-	1.0	2
chr1	11343988	11343989	.	.	-	0.0	2
chr1	11344168	11344169	.	.	-	1.0	2
chr1	11344219	11344220	.	.	-	1.0	2
chr1	11344287	11344288	.	.	-	0.0	2
```

## Command options
```angular2html
tss_eval.py -h

usage: tss_eval (NANOME) [-h] [-v] --dsname DSNAME --runid RUNID --calls CALLS
                         [CALLS ...] [--bgtruth BGTRUTH] [--sep SEP]
                         [--processors PROCESSORS] [-o O] [--enable-cache]
                         [--using-cache] [--output-unified-format]
                         [--chrs CHRS [CHRS ...]] [--tagname TAGNAME]
                         [--debug]

Export read/site level methylation results of all nanopore tools in nanome
paper

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --dsname DSNAME       dataset name
  --runid RUNID         running prefix
  --calls CALLS [CALLS ...]
                        all ONT call results <tool-name>:<file-name> seperated
                        by spaces
  --bgtruth BGTRUTH     background truth file <encode-type>:<file-
                        name1>;<file-name2>
  --sep SEP             seperator for output csv file
  --processors PROCESSORS
                        running processors
  -o O                  output dir
  --enable-cache        if enable cache functions
  --using-cache         if use cache files
  --output-unified-format
                        True:output read level results(1-based start), False:
                        output site-level results(0-based start)
  --chrs CHRS [CHRS ...]
                        chromosome list
  --tagname TAGNAME     output unified file tagname
  --debug               if output debug info
```
