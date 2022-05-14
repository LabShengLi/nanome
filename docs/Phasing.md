**This is an explanation of how to use NANOME pipeline phasing and 5hmC call. For general usage, please check [Usage](https://github.com/LabShengLi/nanome/blob/master/docs/Usage.md).**

## Phasing reads for methylaiton
NANOME support phasing long-reads into haplotypes using [Clair3](https://github.com/HKU-BAL/Clair3) (variant calling) and [Whatshap](https://whatshap.readthedocs.io/en/latest/) (long-read phasing). The option for running phasing step is `--phasing`.

Below is an example:
```angular2html
nextflow run LabShengLi/nanome\
    -profile singularity,winter \
    --dsname NA19240_PhasingDemo \
    --input https://storage.googleapis.com/jax-nanopore-01-project-data/nanome-input/NA19240_phasing_test_data_chr11_chr15_n49.tar.gz\
    --genome hg38 --phasing --hmc --ctg_name 'chr11,chr15'
```

Phasing outputs is below, methylation outputs for Megalodon and NANOME will be split into H1 and H2 haplotypes:
```angular2html
[dd/10cf26] process > EnvCheck (NA19240_PhasingDemo)      [100%] 1 of 1 ✔
[5c/6522a7] process > Untar (NA19240.tar)     [100%] 1 of 1 ✔
[f1/e53635] process > Basecall (NA19240.tar)  [100%] 1 of 1 ✔
[60/73dd05] process > QCExport (NA19240_PhasingDemo)      [100%] 1 of 1 ✔
[41/e623cd] process > Megalodon (NA19240.tar) [100%] 1 of 1 ✔
[c9/65ae15] process > MgldnComb (NA19240_PhasingDemo)     [100%] 1 of 1 ✔
[18/91c19b] process > Report (NA19240_PhasingDemo)        [100%] 1 of 1 ✔
[06/93d747] process > Clair3 (NA19240_PhasingDemo)        [100%] 1 of 1 ✔
[ff/8f16be] process > Phasing (NA19240_PhasingDemo)       [100%] 1 of 1 ✔
Completed at: 26-Apr-2022 22:28:00
Duration    : 2m 38s
CPU hours   : 0.5 (68.3% cached)
Succeeded   : 3
Cached      : 6

ls results/NA19240_PhasingDemo-phasing/
```

## Call 5hmC by Megalodon
NANOME support 5hmC call by Megalodon using param `--hmc`.
