suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(R.utils))

data.table::setDTthreads(4L) 

column_names <-c('Chrom', 'Start', 'End', 'Strand', 'Dataset', 'BSseq_freq', 'BSseq_cov', 
                 'Nanopolish_freq', 'Nanopolish_cov','Megalodon_freq', 'Megalodon_cov','DeepSignal_freq','DeepSignal_cov', 
                 'Guppy_freq','Guppy_cov','Tombo_freq', 'Tombo_cov','DeepMod_freq', 'DeepMod_cov','METEORE_freq', 'METEORE_cov',
                 'Singleton', 'Genomic_location', 'CpG_location','CG_density','Repetitive_regions')
df <- data.table::fread(input='./data/total.unique.tsv.bz2', sep="\t")

data.table::setnames(df, column_names)
df

chr_order <- factor(df$Chrom, levels = c(paste("chr",1:22,sep=""),"chrX","chrY"))
dataset_order <- factor(df$Dataset, levels = c("NA19240", "NA12878", "APL","K562"))
strand_order <- factor(df$Strand, levels = c("+","-"))
singleton_order <- factor(df$Singleton, levels = c("singleton","non-singleton"))
genomic_order <- factor(df$Genomic_location, levels = c("promoter", "exon", "intron", "intergenic"))
CGdensity_order <- factor(df$CG_density, levels = c("100%", "80%",  "60%", "40%", "20%"))

df %>% dplyr::mutate_if(is.character,funs(factor(.)))

options(DT.options = list(columnDefs = list(list(className = 'dt-center', width = '6%', targets = "all")),
                          pageLength = 25,  
                          lengthMenu = c(25, 50, 75),
                          scrollX = T))

