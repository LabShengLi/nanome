**This is an explanation of how to use NANOME pipeline phasing and 5hmC call. For general usage, please check [Usage](https://github.com/TheJacksonLaboratory/nanome/blob/master/docs/Usage.md).**

## Phasing reads for methylaiton
NANOME support phasing long-reads into H1 and H2 using [Clair3](https://github.com/HKU-BAL/Clair3) (variant calling) and [Whatshap](https://whatshap.readthedocs.io/en/latest/) (long-read phasing). The option for running phasing step is `--phasing`.

Below is an example:
```angular2html
nextflow run TheJacksonLaboratory/nanome\
    -profile singularity,winter \
    --dsname NA12878_CHR22_200 \
    --input https://storage.googleapis.com/jax-nanopore-01-project-data/nanome-input/na12878_chr22_200.tar.gz\
    --genome hg38_chr22\
    --ctg_name chr22\
    --runGuppy false --runNanopolish false --runDeepSignal false\
    --phasing
```

Phasing outputs is below, methylation outputs for Megalodon and NANOME will be split into H1 and H2 haplotypes:
```angular2html
[dd/10cf26] process > EnvCheck (NA12878_CHR22_200)      [100%] 1 of 1 ✔
[5c/6522a7] process > Untar (na12878_chr22_200.tar)     [100%] 1 of 1 ✔
[f1/e53635] process > Basecall (na12878_chr22_200.tar)  [100%] 1 of 1 ✔
[60/73dd05] process > QCExport (NA12878_CHR22_200)      [100%] 1 of 1 ✔
[41/e623cd] process > Megalodon (na12878_chr22_200.tar) [100%] 1 of 1 ✔
[c9/65ae15] process > MgldnComb (NA12878_CHR22_200)     [100%] 1 of 1 ✔
[18/91c19b] process > Report (NA12878_CHR22_200)        [100%] 1 of 1 ✔
[06/93d747] process > Clair3 (NA12878_CHR22_200)        [100%] 1 of 1 ✔
[ff/8f16be] process > Phasing (NA12878_CHR22_200)       [100%] 1 of 1 ✔
Completed at: 26-Apr-2022 22:28:00
Duration    : 2m 38s
CPU hours   : 0.5 (68.3% cached)
Succeeded   : 3
Cached      : 6

ls results/NA12878_CHR22_200-phasing/
NA12878_CHR22_200_clair3_out  NA12878_CHR22_200_hp_split

ls results/NA12878_CHR22_200-phasing/NA12878_CHR22_200_hp_split
NA12878_CHR22_200_megalodon_read_pred_chr22_H1.tsv.gz
NA12878_CHR22_200_megalodon_read_pred_chr22_H2.tsv.gz
NA12878_CHR22_200_megalodon_site_freq_chr22_H1.tsv.gz
NA12878_CHR22_200_megalodon_site_freq_chr22_H2.tsv.gz
NA12878_CHR22_200_nanome_na12878_xgboostna3t_read_pred_chr22_H1.tsv.gz
NA12878_CHR22_200_nanome_na12878_xgboostna3t_read_pred_chr22_H2.tsv.gz
NA12878_CHR22_200_nanome_na12878_xgboostna3t_site_freq_chr22_H1.tsv.gz
NA12878_CHR22_200_nanome_na12878_xgboostna3t_site_freq_chr22_H2.tsv.gz
```

## Call 5hmC by Megalodon
NANOME support 5hmC call by Megalodon using param `--hmc`.
