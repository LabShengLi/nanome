## NANOME pipeline options


| Name | Default                                   | Description                                                 |
|------|-------------------------------------------|-------------------------------------------------------------|
| dsname   | NA                                        | Analyzed dataset name, **required**                             |
| input   | NA                                        | Input fast5 files directory or compressed tar/tar.gz file, **required**   |
| genome   | hg38                                      | Reference genome                                            |
| outdir   | results                                   | Output directory                                            |
| cleanup   | false                                     | Clean Nextflow work/ directory if true                      |
| runNanopolish   | true                                      | If running Nanpolish                                        |
| runMegalodon   | true                                      | If running Megalodon                                        |
| runDeepSignal   | true                                      | If running DeepSignal                                       |
| runGuppy   | false                                     | If running Guppy methylation-calling                        |
| runNANOME   | true                                      | If running NANOME consensus model                           |
| hmc   | false                                     | If use 5hmC_5mC model for Megalodon                         |
| phasing   | false                                     | If perform methylation phasing                              |
| ctg_name   | null                                      | Regions for variant calling and phasing, such as chr1, etc. |
| GUPPY_BASECALL_MODEL   | dna_r9.4.1_450bps_hac.cfg                 | Guppy basecalling model                                     |
| GUPPY_METHCALL_MODEL   | dna_r9.4.1_450bps_modbases_5mc_cg_hac.cfg | Guppy methylation-calling model                             |
| remoraModel   | dna_r9.4.1_e8                             | Megalodon model                                             |
| GUPPY_TIMEOUT   | 500                                       | Guppy timeout parameter for Megalodon                       |
| NANOME_MODEL   | NANOME3T                                  | NANOME consensus model, 'NANOME3T' or 'NANOME2T'            |
| DEEPSIGNAL_MODEL   | bn_17.sn_360.epoch_9.ckpt                 | DeepSignal model                                            |

