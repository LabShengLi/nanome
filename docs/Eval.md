**This is an explanation of how to use script to perform read-level and site-level performance evaluation.**

**Please note that we strongly suggest providing the whole genome-wide CpGs for performance comparison across tools, we do not suggest selecting a portion of CpGs for comparison**.

User needs to provide two kind of files as input for our proposed standardized benchmarking framework:
1. ONT tools' methylation-calling raw outputs,
1. BS-seq data. 
 
The [genome-annotation](https://storage.googleapis.com/jax-nanopore-01-project-data/nanome_paper/genome-annotation.tar.gz) files may be needed if performances at specific genomic regions are interested. 

You can install NANOME standardized genome-wide evaluation tool (**nanome-jax**) from [PyPI](https://pypi.org/project/nanome-jax/):
```angular2html
# Creat conda environment and install required packages
conda create --name py39 python=3.9
conda activate py39
conda install -c anaconda scikit-learn cython
conda install -c conda-forge -c bioconda pybedtools

# Install nanome-jax for genome-wide evaluation
pip install nanome-jax
pip show nanome-jax
```

**If your system support docker or singularity, you can directly run evaluation tool without any installations**:
```angular2html
## Check nanome-jax package in docker
docker run liuyangzzu/nanome pip show nanome-jax

## Check nanome-jax package in singularity 
singularity exec docker://liuyangzzu/nanome pip show nanome-jax
```

# 1. Read-level performance evaluation on common/shared CpGs
## Script for read-level evaluation
The script `read_level_eval.py` is designed for general purpose of read-level performance evaluation for all kinds of tools on common CpGs reported by tools and BS-seq.

```angular2html
read_level_eval.py -v
```

## Sample usage
User needs to specify methylation-calling raw results of each tool and BS-seq replicates to execute the read-level evaluation script. 

`--calls` params are raw input of all tools, format is like ` <tool-name>:<file-encode>:<file-name>`, where `tool-name` can be any name for your raw input, `file-encode` is the input format, can be Nanopolish, Megalodon, DeepSignal, Guppy, Tombo, METEORE, DeepMod, NANOME, and `file-name` is the location of raw input. `--bgtruth` is the background truth param like `<encode-type>:<file-name1>;<file-name2>`, where `encode-type` can be 'encode' or 'bismark', and we support at most two replicates file inputs now.
```angular2html
read_level_eval.py \
    --dsname TestData\
    --runid MethPerf-TestData_RRBS_2Reps\
    --calls\
        Nanopolish:Nanopolish:[Nanopolish-calls]\
        Megalodon:Megalodon:[Megalodon-calls]\
        DeepSignal:DeepSignal:[DeepSignal-calls]\
    --bgtruth [encode-type]:[BS-seq-rep1];[BS-seq-rep2]\
    --genome-annotation [genome-annotation]
```
Sample results can be found at [read-level outputs](https://github.com/LabShengLi/nanome/blob/master/docs/resources/read_level_output.txt).

For using docker or singularity, the command prefix is like below:
```angular2html
# Docker style usage
docker run -v $PWD:$PWD -w $PWD -it liuyangzzu/nanome read_level_eval.py

# Singularity style usage
singularity exec -e docker://liuyangzzu/nanome read_level_eval.py
```

## Command options
```angular2html
read_level_eval.py -h

usage: read_level_eval (NANOME) [-h] [-v] --dsname DSNAME --runid RUNID
                                --calls CALLS [CALLS ...] [--bgtruth BGTRUTH]
                                [--genome-annotation GENOME_ANNOTATION]
                                [--min-bgtruth-cov MIN_BGTRUTH_COV]
                                [--toolcov-cutoff TOOLCOV_CUTOFF]
                                [--processors PROCESSORS] [--report-no-join]
                                [--chrSet CHRSET [CHRSET ...]] [-o O]
                                [--enable-cache] [--using-cache]
                                [--distribution] [--bsseq-report]
                                [--analysis ANALYSIS] [--save-curve-data]
                                [--large-mem] [--bedtools-tmp BEDTOOLS_TMP]
                                [--cache-dir CACHE_DIR] [--disable-bed-check]
                                [--mpi] [--mpi-import] [--config] [--verbose]

Read-level performance evaluation in nanome paper

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --dsname DSNAME       dataset name
  --runid RUNID         running prefix/output folder name, such as MethPerf-
                        Dataset_WGBS_2Reps
  --calls CALLS [CALLS ...]
                        all ONT call results <tool-name>:<file-encode>:<file-
                        name>[<cutoff0>:<cutoff1>] seperated by space, tool-
                        name/file-encode can be Nanopolish, Megalodon,
                        DeepSignal, Guppy, Tombo, METEORE, DeepMod, NANOME.
                        The optional cutoff0 and cutoff1 can be specified if
                        user wants to change default cutoff values, i.e.,
                        Nanpolish uses -2:2, Megalodon uses 0.2:0.8, etc.
  --bgtruth BGTRUTH     background truth file <encode-type>:<file-
                        name1>;<file-name2>, encode-type can be 'encode' or
                        'bismark'
  --genome-annotation GENOME_ANNOTATION
                        genome annotation dir, contain BED files
  --min-bgtruth-cov MIN_BGTRUTH_COV
                        min bg-truth coverage cutoff, default is 5
  --toolcov-cutoff TOOLCOV_CUTOFF
                        cutoff for coverage in nanopore tools, default is >=1
  --processors PROCESSORS
                        number of processors used, default is 1
  --report-no-join      true if report not on joined sets
  --chrSet CHRSET [CHRSET ...]
                        chromosome list, default is human chr1-22, X and Y
  -o O                  output base dir
  --enable-cache        if enable cache functions
  --using-cache         if use cache files
  --distribution        if report singleton/nonsingleton distributions at all
                        regions
  --bsseq-report        if report singleton/nonsingleton in bs-seq
  --analysis ANALYSIS   special analysis specifications for ecoli
  --save-curve-data     if save pred/truth points for curve plot
  --large-mem           if using large memory (>100GB) for speed up
  --bedtools-tmp BEDTOOLS_TMP
                        bedtools temp dir
  --cache-dir CACHE_DIR
                        cache dir used for loading calls/bs-seq (speed up
                        running)
  --disable-bed-check   if disable auto-checking the 0/1 base format for
                        genome annotations
  --mpi                 if using multi-processing/threading for evaluation, it
                        can speed-up but may need more memory
  --mpi-import          if using multi-processing/threading for import, it can
                        speed-up, only for small size data
  --config              if print out config file for genome annotation
  --verbose             if output verbose info
```

# 2. Site-level performance evaluation on common/shared CpGs

## Script for site-level evaluation
The script `site_level_eval.py` is designed for general purpose of site-level performance evaluation for all kinds of tools on common CpGs reported by tools and BS-seq.

```angular2html
site_level_eval.py -v
```


## Sample usage
User needs to specify methylation-calling raw results of each tool and BS-seq replicates to execute the site-level evaluation script.

```angular2html
site_level_eval.py \
    --dsname TestData\    
    --runid MethCorr-TestData_RRBS_2Reps\
    --calls\
        Nanopolish:Nanopolish:[Nanopolish-calls]\
        Megalodon:Megalodon:[Megalodon-calls]\
        DeepSignal:DeepSignal:[DeepSignal-calls]\
    --bgtruth [encode-type]:[BS-seq-rep1];[BS-seq-rep2]\
    --genome-annotation [genome-annotation]\
    --min-bgtruth-cov 3 --toolcov-cutoff 1\
    --beddir [read-level-analysis-dir]
```

Sample results can be found at [site-level outputs](https://github.com/LabShengLi/nanome/blob/master/docs/resources/site_level_output.txt).

## Command options
```angular2html
site_level_eval.py -h

usage: site_level_eval (NANOME) [-h] [-v] --dsname DSNAME --runid RUNID --calls CALLS [CALLS ...] [--bgtruth BGTRUTH]
                                [--genome-annotation GENOME_ANNOTATION] [--beddir BEDDIR]
                                [--min-bgtruth-cov MIN_BGTRUTH_COV] [--toolcov-cutoff TOOLCOV_CUTOFF]
                                [--chrSet CHRSET [CHRSET ...]] [--sep SEP] [--processors PROCESSORS] [-o O] [--gen-venn]
                                [--sort] [--deduplicate] [--summary-coverage] [--region-coe-report] [--report-no-join]
                                [--enable-cache] [--using-cache] [--plot] [--bedtools-tmp BEDTOOLS_TMP]
                                [--cache-dir CACHE_DIR] [--large-mem] [--disable-bed-check] [--mpi] [--config] [--verbose]

Site-level correlation analysis in nanome paper

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --dsname DSNAME       dataset name
  --runid RUNID         running prefix/output folder name, such as MethCorr-Dataset_WGBS_2Reps
  --calls CALLS [CALLS ...]
                        all ONT call results <tool-name>:<file-encode>:<file-name> seperated by spaces, tool-name/file-
                        encode can be Nanopolish, Megalodon, DeepSignal, Guppy, Tombo, METEORE, DeepMod, NANOME. The
                        optional cutoff0 and cutoff1 can be specified if user wants to change default cutoff values, i.e.,
                        Nanpolish uses -2:2, Megalodon uses 0.2:0.8, etc.
  --bgtruth BGTRUTH     background truth file <encode-type>:<file-name1>;<file-name2>, encode-type can be 'encode' or
                        'bismark'
  --genome-annotation GENOME_ANNOTATION
                        genome annotation dir, contain BED files
  --beddir BEDDIR       base dir for concordant/discordant BED files generated by read-level analysis, make sure the dsname
                        is same
  --min-bgtruth-cov MIN_BGTRUTH_COV
                        cutoff for coverage in bg-truth, default is >=5
  --toolcov-cutoff TOOLCOV_CUTOFF
                        cutoff for coverage in nanopore tools, default is >=3
  --chrSet CHRSET [CHRSET ...]
                        chromosome list, default is human chr1-22, X and Y
  --sep SEP             seperator for output csv file
  --processors PROCESSORS
                        number of processors used, default is 1
  -o O                  output base dir
  --gen-venn            if generate CpGs files for venn data analysis
  --sort                if sort outputs
  --deduplicate         if deduplicate not unique records needed
  --summary-coverage    if summarize coverage at each region
  --region-coe-report   if report PCC value at each region
  --report-no-join      if output no-join report also
  --enable-cache        if enable cache functions
  --using-cache         if use cache files
  --plot                if plot the correlation matrix figure
  --bedtools-tmp BEDTOOLS_TMP
                        bedtools temp dir
  --cache-dir CACHE_DIR
                        cache dir used for loading calls/bs-seq(speed up running)
  --large-mem           if using large memory (>100GB) for speed up
  --disable-bed-check   if disable auto-checking the 0/1 base format for genome annotations
  --mpi                 if using multi-processing/threading for evaluation, it can speed-up but need more memory
  --config              if print out config file for genome annotation
  --verbose             if output verbose info
```

# 3. Unified read/site level format conversion


## Script for unified format conversion
The script `tss_eval.py` is designed for general purpose of converting raw results of all kinds of tools into a unified format.

```angular2html
tss_eval.py -v
```

## Sample usage for read-level format unification
The option `--read-level-format` is used for generate unified read-level format for ONT tools.
```angular2html
tss_eval.py \
    --dsname [dataset-name]\
    --runid Read_Level-[dataset-name] \
    --calls \
        [tool-name]:[encode-name]:[tool-call-filename] \
    --read-level-format\
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
The columns for read level output are:
1. read-id, 
2. chromosome, 
3. position (1-based), 
4. strand, 
5. score of log-ratio for probability of 5mC vs. 5C.

## Sample usage for site-level format unification
Below will generate the site level unified format for tools/BS-seq data.
```angular2html
tss_eval.py \
    --dsname [dataset-name]\
    --runid Site_Level-[dataset-name]\
    --bgtruth [encode-type]:[BS-seq-file1];[BS-seq-file2]\
    --calls \
        [tool-name]:[encode-type]:[tool-call-filename]\
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

The columns for site level output are:
1. chromosome, 
2. start (0-based), 
3. end (1-based), 
4. NA, 
5. NA, 
6. strand, 
7. methylation frequency, 
8. coverage.

## Command options
```angular2html
tss_eval.py -h

usage: tss_eval (NANOME) [-h] [-v] --dsname DSNAME --runid RUNID [--calls CALLS [CALLS ...]] [--bgtruth BGTRUTH]
                         [--read-level-format] [--sep SEP] [--processors PROCESSORS] [-o O] [--enable-cache] [--using-cache]
                         [--sort] [--deduplicate] [--cache-dir CACHE_DIR] [--chrSet CHRSET [CHRSET ...]] [--tagname TAGNAME]
                         [--verbose]

Export read/site level methylation results of all nanopore tools in nanome paper

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --dsname DSNAME       dataset name
  --runid RUNID         running prefix/output dir name
  --calls CALLS [CALLS ...]
                        all ONT call results <tool-name>:<file-encode>:<file-name>[<cutoff0>:<cutoff1>] seperated by space,
                        tool-name/file-encode can be Nanopolish, Megalodon, DeepSignal, Guppy, Tombo, METEORE, DeepMod,
                        NANOME. The optional cutoff0 and cutoff1 can be specified if user wants to change default cutoff
                        values, i.e., Nanpolish uses -2:2, Megalodon uses 0.2:0.8, etc.
  --bgtruth BGTRUTH     background truth file <encode-type>:<file-name1>;<file-name2>
  --read-level-format   if true, it will output read level results (1-based start), else it will output site-level results
                        (0-based start, 1-based end)
  --sep SEP             seperator for output csv file, default is tab character
  --processors PROCESSORS
                        running processors, default is 1
  -o O                  output base dir
  --enable-cache        if enable cache functions
  --using-cache         if use cache files
  --sort                if sort bed output
  --deduplicate         if deduplicate not unique records
  --cache-dir CACHE_DIR
                        cache dir used for loading calls/bs-seq (speed up running)
  --chrSet CHRSET [CHRSET ...]
                        chromosome list, default is human chromosome chr1-22, X and Y
  --tagname TAGNAME     output unified file's tagname
  --verbose             if output verbose info
```

# 4. Intersect your data with genomic regions

For specific genomic region analysis, below is the tool to intersect data file with specific region BED file.

```angular2html
region_intersect.py -h

usage: region_intersect (NANOME) [-h] [-v] -i I -r R -o O [--sep SEP] [--bed-format BED_FORMAT] [--strand-intersect]
                                 [--header HEADER] [--verbose]

Intersect data with region file

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -i I                  input data file
  -r R                  region BED file
  -o O                  output file
  --sep SEP             file seperator, default is TAB
  --bed-format BED_FORMAT
                        BED file start format 0/1, default is 1
  --strand-intersect    if BED file intersect using strand sensitive mode
  --header HEADER       if input file contain header, set to 0, default is None
  --verbose             if output verbose info
```

We provide a list of regions on GCP storage: https://storage.googleapis.com/jax-nanopore-01-project-data/nanome_paper/genome-annotation.tar.gz .
