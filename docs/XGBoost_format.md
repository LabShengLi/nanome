**This is the input/output format for XGBoost model train and prediction.**

## Train input/output format
Below are the example input, columns are:
1. **Chr**  
   Chromosome where the CpG site is located (e.g., chr1, chrX).

2. **Pos**  
   Genomic coordinate (1-based) of the CpG site.

3. **Strand**  
   DNA strand of the CpG site, either '+' (forward) or '−' (reverse).

4. **Pos_deprecated**  
   Legacy or alternative position (may refer to older genome assemblies or internal indexing).

5. **ID**  
   Read ID.

6. **read_strand**  
   Strand orientation.

7. **k_mer**  
   Sequence context surrounding the CpG site.

8. **methy_label**  
   Ground truth label (0 = unmethylated, 1 = methylated).

9. **megalodon**  
   Methylation prediction score or binary call from the **Megalodon** tool.

10. **nanopolish**  
    Methylation score or binary call from **Nanopolish**.

11. **deepsignal**  
    Methylation prediction from the **DeepSignal** model.

12. **Freq**  
    Frequency of methylation at this site across all reads (methylated reads / total reads) in BS-seq.

13. **Coverage**  
    Total number of reads covering this CpG site in BS-seq.

14. **Label**  
    Binary or categorical label for evaluation.

15. **Region**  
    Genomic annotation for the site (e.g., singleton, concordant, etc.).


```angular2html
Chr     Pos     Strand  Pos_deprecated  ID      read_strand     k_mer  methy_label     megalodon     nanopolish      deepsignal     Freq    Coverage        Label   Region
chr1    11002684        -       -1      790a97af-c74d-44f0-bb19-b06e1e90eaf9    t       ACCACACCCGGCTAATT     1      -1.77   -0.1073302887781659     5.972520033294956       0.0     23      0       Singleton
```

Output is the trained model (`.pkl` file) and cross validation report (`.xlsx` file), example output files are:

```angular2html
NA12878_chr1_xgboost_basic_w_seq_cv_results.xlsx
NA12878_chr1_xgboost_basic_w_seq_model.pkl
```


## Predict script input/output format

Predict script contains input files from three tools' read-level prediction results, each prediction result file format is below:

```angular2html
ID      Chr     Pos     Strand  Score
48b198e0-1733-40d1-b23e-60b82033e3be    chr1    9037301 +       -1.3768708366148037
48b198e0-1733-40d1-b23e-60b82033e3be    chr1    9037318 +       -1.5710167549897966
```

Columns are:

1. **ID**  
   Read ID.

2. **Chr**  
   Chromosome (e.g., chr1, chrX).

3. **Pos**  
   Genomic position (1-based location).

4. **Strand**  
   DNA strand, either '+' (forward) or '−' (reverse).

For train model with DNA sequence features, the feature input is also needed for below format (same as DeepSignal feature extraction format):

```angular2html
chr17   3461708 +       -1      fabc96c5-2543-4c08-82fc-b4ca7b5023c7    t       AAAAAATACGATATCTT     0.099924,-0.204013,-0.38477,-0.574565,-0.961772,-0.824376,1.259275,0.26587,0.29561,-1.13664,-1.124149,1.141637,-0.162377,1.992243,0.312264,0.387207,0.774414  0.103503,0.119515,0.093031,0.113566,0.082928,0.117615,0.33196,0.205535,0.22629,0.086236,0.046735,0.224357,0.556972,0.118441,0.549111,0.131809,0.065699      3,3,41,3,25,3,11,7,3,6,3,5,5,6,3,25,6   0.0,0.0,0.0,0.0,0.0,0.0,0.23732,0.074943,-0.012491,0.0,0.0,0.0,0.0,0.0,0.0,0.0;0.0,0.0,0.0,0.0,0.0,0.0,-0.037472,-0.262302,-0.312264,0.0,0.0,0.0,0.0,0.0,0.0,0.0;-0.362226,-0.399698,-0.349735,-0.324754,-0.374717,-0.337245,-0.387207,-0.44966,-0.299773,-0.549584,-0.387207,-0.324754,-0.299773,-0.412188,-0.474641,-0.46215;0.0,0.0,0.0,0.0,0.0,0.0,-0.44966,-0.549584,-0.724452,0.0,0.0,0.0,0.0,0.0,0.0,0.0;-0.89932,-0.961772,-1.061697,-1.049206,-0.999244,-0.949282,-0.824376,-0.961772,-0.849357,-1.036716,-0.861848,-1.024225,-1.061697,-0.949282,-0.986754,-0.961772;0.0,0.0,0.0,0.0,0.0,0.0,-0.986754,-0.711961,-0.774414,0.0,0.0,0.0,0.0,0.0,0.0,0.0;0.0,0.0,0.299773,1.486376,1.511357,1.34898,1.511357,1.299017,1.436413,1.074187,1.336489,1.149131,1.398942,0.0,0.0,0.0;0.0,0.0,0.0,0.0,0.699471,0.287283,0.187358,0.249811,-0.012491,0.124906,0.324754,0.0,0.0,0.0,0.0,0.0;0.0,0.0,0.0,0.0,0.0,0.0,0.549584,0.337245,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0;0.0,0.0,0.0,0.0,0.0,-1.13664,-1.036716,-1.049206,-1.299017,-1.13664,-1.161621,0.0,0.0,0.0,0.0,0.0;0.0,0.0,0.0,0.0,0.0,0.0,-1.074187,-1.111659,-1.186602,0.0,0.0,0.0,0.0,0.0,0.0,0.0;0.0,0.0,0.0,0.0,0.0,0.786905,1.061697,1.473885,1.149131,1.236565,0.0,0.0,0.0,0.0,0.0,0.0;0.0,0.0,0.0,0.0,0.0,0.811886,-0.524603,-0.574565,-0.637018,0.112415,0.0,0.0,0.0,0.0,0.0,0.0;0.0,0.0,0.0,0.0,0.0,1.998488,1.748677,2.04845,2.135884,2.023469,1.998488,0.0,0.0,0.0,0.0,0.0;0.0,0.0,0.0,0.0,0.0,0.0,1.024225,-0.312264,0.22483,0.0,0.0,0.0,0.0,0.0,0.0,0.0;0.324754,0.599546,0.324754,0.287283,0.162377,0.46215,0.262302,0.487131,0.374717,0.412188,0.549584,0.349735,0.524603,0.324754,0.637018,0.512113;0.0,0.0,0.0,0.0,0.0,0.649509,0.849357,0.786905,0.799395,0.824376,0.736943,0.0,0.0,0.0,0.0,0.0   1
```

Each column is below:
1. **Chr**  
   Chromosome name (e.g., chr17) where the methylation site is located.

2. **Pos**  
   Genomic coordinate of the cytosine site (typically the center base in the k-mer context).

3. **Strand**  
   DNA strand orientation of the methylation site: '+' (forward) or '−' (reverse).

4. **Pos_deprecated**  
   Deprecated or placeholder position; often unused or set to -1 in newer DeepSignal outputs.

5. **ID**  
   Read ID.

6. **read_strand**  
   Strand of the sequencing read (usually 't' for template or 'c' for complement).

7. **k_mer**  
   DNA sequence context (typically 17-mer) surrounding the methylation site, centered at the cytosine.

8. **signal_means**  
   Comma-separated list of normalized mean signal values for each base in the k-mer context.

9. **signal_stds**  
   Comma-separated list of normalized standard deviations of raw signals for each base in the k-mer context.

10. **signal_lengths**  
    Comma-separated list of the number of raw signal samples aligned to each base in the k-mer.

11. **event_features**  
    A long semicolon-separated matrix-like field where each segment contains comma-separated DeepSignal event-level features (e.g., raw signal summaries, deltas, etc.) for each base in the k-mer.  
    These features are extracted from the raw nanopore signal for DeepSignal’s deep learning model. Each row corresponds to a different feature type (e.g., dwell time, delta current).

12. **methy_label**  
    The ground truth label for supervised learning or benchmarking:  
    - `1` = methylated  
    - `0` = unmethylated  
    - Sometimes `-1` indicates unknown or not labeled.


Predict script output file format is below:
```angular2html
ID      Chr     Pos     Strand  megalodon       nanopolish      deepsignal      Prediction   Prob_methylation
13766665-3638-4223-931e-733f3e6e9242    chr1    3605228 +       -3.0417192588801085     -2.51-1.204658646447991       0       0.08439023888455542
```

1. **ID**  
    Read ID.
2. **Chr**  
   Chromosome where the CpG site is located (e.g., chr1, chr17).

3. **Pos**  
   Genomic coordinate (1-based position) of the cytosine site.

4. **Strand**  
   DNA strand of the CpG site: '+' (forward) or '−' (reverse).

5. **megalodon**  
   Methylation prediction output.

6. **nanopolish**  
   Methylation prediction from **Nanopolish**.

7. **deepsignal**  
   Methylation prediction from **DeepSignal**.

8. **Prediction**  
   Final consensus or model-inferred binary classification:  
   - `1` = methylated  
   - `0` = unmethylated

9. **Prob_methylation**  
   Final consensus or model-inferred probability that the site is methylated, typically between 0 and 1.
