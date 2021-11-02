suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(R.utils))
data.table::setDTthreads(4L) 

column_names <-c('Chrom', 'Start', 'End', 'Strand','Dataset', 'BSseq_freq', 'BSseq_cov', 
                 'Nanopolish_freq', 'Nanopolish_cov','Megalodon_freq', 'Megalodon_cov','DeepSignal_freq','DeepSignal_cov', 
                 'Guppy_freq','Guppy_cov','Tombo_freq', 'Tombo_cov','DeepMod_freq', 'DeepMod_cov','METEORE_freq', 'METEORE_cov',
                 'Singleton', 'Genomic_location', 'CpG_location','CGdensity','Repetitive_region')

df <- data.table::fread(input='./data/total.tsv.bz2', sep="\t")
data.table::setnames(df, column_names)
df

#unique(df$Genome_location)
#https://stackoverflow.com/questions/47539860/r-remove-duplicates-from-a-dataframe-based-on-categories-in-a-column
GenomeLocationOrder <-  c("promoter", "exon", "intron", "intergenic")
df_new <- df
df_new$Genomic_location <- factor(df_new$Genomic_location, levels=GenomeLocationOrder)
df_new <- df_new[order(df_new$Genomic_location),]
# remove multiple genomic features
df_new <- df_new[!duplicated(df_new[,c('Chrom', 'Start', 'End', 'Strand','Dataset', 'BSseq_freq', 'BSseq_cov', 
                                       'Nanopolish_freq', 'Nanopolish_cov','Megalodon_freq', 'Megalodon_cov','DeepSignal_freq','DeepSignal_cov', 
                                       'Guppy_freq','Guppy_cov','Tombo_freq', 'Tombo_cov','DeepMod_freq', 'DeepMod_cov','METEORE_freq', 'METEORE_cov',
                                       'Singleton', 'Genomic_location', 'CpG_location','CGdensity','Repetitive_region')]),]

dataset_order <- c("NA19240", "NA12878", "APL", "K562")
df_new %>%
  arrange(factor(Dataset, levels = dataset_order),
          factor(Chrom, levels = c(paste("chr",1:22,sep=""),"chrX","chrY")), 
          Start)
dim(df)
dim(df_new)
df_new
summary(as.factor(df_new$Dataset))
summary(as.factor(df_new$Singleton))
summary(as.factor(df_new$Genomic_location))
summary(as.factor(df_new$CpG_location))
summary(as.factor(df_new$CGdensity))
summary(as.factor(df_new$Repetitive_region))

write.table(df_new, file='./data/total.unique.tsv', sep="\t", row.names = FALSE, col.names=FALSE, quote = FALSE)

#Check table
column_names <-c('Chrom', 'Start', 'End', 'Strand', 'Dataset', 'BSseq_freq', 'BSseq_cov', 
                 'Nanopolish_freq', 'Nanopolish_cov','Megalodon_freq', 'Megalodon_cov','DeepSignal_freq','DeepSignal_cov', 
                 'Guppy_freq','Guppy_cov','Tombo_freq', 'Tombo_cov','DeepMod_freq', 'DeepMod_cov','METEORE_freq', 'METEORE_cov',
                 'Singleton', 'Genomic_location', 'CpG_location','CG_density','Repetitive_regions')
df <- data.table::fread(input='./data/total.unique.tsv.bz2', sep="\t")

summary(as.factor(df$Singleton))
summary(as.factor(df$Genomic_location))
summary(as.factor(df$CpG_location))
summary(as.factor(df$CG_density))
summary(as.factor(df$Repetitive_region))
summary(as.factor(df$Dataset))
summary(as.factor(df$Strand))
summary(as.factor(df$Chrom))


<b>Citation</b><br>
  If you used the dataset for your research, please cite the following publication:
  <br>
  Liu, Y., Rosikiewicz, W., Pan, Z. et al. 
DNA methylation-calling tools for Oxford Nanopore sequencing: a survey and human epigenome-wide evaluation. 
Genome Biol 22, 295 (2021).,
<a href=paperlink title="DOI" target="_blank">
  https://doi.org/10.1186/s13059-021-02510-z</a><br>
  <br>

  <br/>
  The data presented here associated with the publication/preprint in the Citation section.