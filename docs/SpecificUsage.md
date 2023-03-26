**This is an explanation of how to use NANOME pipeline on some specific scenarios. For general usage, please check [Usage](https://github.com/LabShengLi/nanome/blob/master/docs/Usage.md).**
## 1. How to use other reference genome in NANOME?

NANOME support other reference genome. Below is an example of running NANOME for any other reference genomes, please make sure put reference genome file .fasta and the indexed file into directory [reference-genome-dir], the `--chrSet` is the chromosomes params for the specific genome. 

```angular2html
nextflow run LabShengLi/nanome\
    -profile singularity \
    --dsname [dataset-name]\
    --input [input-file]\
    --genome [reference-genome-dir]\
    --chrSet '[chromosomes separated by a space]'
```

The reference genome folder `other_ref` will be like below:
```angular2html
tree other_ref/

other_ref
├── Ecoli_k12_mg1655.fasta
├── Ecoli_k12_mg1655.fasta.fai
```

Note: 
* NANOME support the default behaviour of each tool's running, if the tool performed on the specified genome that is same as human/E. coli data. 
* please ensure there is only one `.fasta` file in genome reference folder, the index file is same name with suffix as `.fasta.fai`
* please set param `--chrSet` with the interested chromosomes names, seperated by a space, such as `--chrSet 'chr1 chr2 chr3'` (the single quote character is important)
* for human or ecoli genome, NANOME will automatically set default chromosomes chr1-22, X, Y for human, NC_000913.3 for ecoli. Users can change the default setting by specifying `--chrSet` param

## 2. How to input basecalled data into NANOME? 

NANOME support the basecalled input data. You can use `--skipBasecall` param to inform NANOME to skip redo the basecalling step. We suggest keeping the Guppy basecalled data/directory structure as input. The `.fastq.gz` file, `sequencing_summary.txt` file, and the `workspace/*.fast5` files are need to be located in the default directory.  

Example:
```angular2html
nextflow run main.nf -profile test,singularity\
    --skipBasecall\
    --input test_data/ecoli_ci_basecalled
```

Basecalled input data directory structure like below:
```angular2html
tree test_data/ecoli_ci_basecalled

test_data/ecoli_ci_basecalled
├── abc.fastq.gz
├── fail
├── pass
├── sequencing_summary.txt
└── workspace
    ├── kelvin_20160617_FN_MN17519_sequencing_run_sample_id_74930_ch138_read698_strand.fast5
    └── kelvin_20160617_FN_MN17519_sequencing_run_sample_id_74930_ch139_read4507_strand.fast5
```
Note:
* please use `.fastq.gz`, not `.fastq`, keep their locations in base folder, or `pass`, `fail` folders
* please keep `sequencing_summary.txt` file for QC analysis
* please keep all fast5 files in `workspace/` directory
* NANOME input can be a directory/tar/tar.gz file, user can compress the basecalled folder using `tar -czf ecoli_ci_basecalled.tar.gz ecoli_ci_basecalled/`
* by default, NANOME run top four performers. If user wants to run Nanopolish tool only, just specify not run other three tools using params `--runMegalodon false --runDeepSignal false --runGuppy false`


## 3. Update latest version of NANOME.
Below is the command for updating the latest version of NANOME, if users have executed NANOME old version previously. `-r [branch-name]` is used to get a specific version from a branch in the GitHub for Nextflow.
```angular2html
nextflow pull LabShengLi/nanome

nextflow pull LabShengLi/nanome -r [branch-name]
```

## 4. Perform only basecall and QC.
Using option `--runMethcall false` will not run methylation calling, it will only perform basecall and QC.

```angular2html
nextflow run LabShengLi/nanome\
    -profile singularity,winter \
    --dsname APL\
    --input '/fastscratch/liuya/nanome/APL_ont_out/APL_sept/sept_dir/*'\
    --genome hg38 \
    --runMethcall false
```

## 5. Support T2T-CHM13 genome

Example as below:
```angular2html
nextflow run $NANOME_DIR \
    -profile test_human,singularity \
    --genome chm13
```

## 6. Support R10.4.1 flow cells

Example as below:
```angular2html
nextflow run $NANOME_DIR \
    -profile test_human,singularity \
    --input https://storage.googleapis.com/jax-nanopore-01-project-data/nanome-input/testdata_r10_4_1.tar.gz \
    --runGuppy \
    --GUPPY_BASECALL_MODEL dna_r10.4.1_e8.2_400bps_hac.cfg \
    --GUPPY_METHCALL_MODEL dna_r10.4.1_e8.2_400bps_modbases_5mc_cg_hac.cfg \
    --runNanopolish false --runDeepSignal false --runMegalodon false
```