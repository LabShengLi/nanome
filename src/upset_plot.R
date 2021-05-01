rm(list = ls())
library(here)
library(UpSetR)
library(eulerr)
library(ComplexHeatmap)
library("optparse")

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "result/venn-data",
              help = "Input dir"),
  make_option(c("-o", "--output"), type = "character", default = getwd(),
              help = "Output dir"),
  make_option(c("--dsname"), type = "character", default = "NA19240",
              help = "Dataset name")
)

opt = parse_args(OptionParser(option_list = option_list))

print(opt)

indir = opt$input
outdir = opt$output
dsname = opt$dsname

data_dir = indir #here('result', 'venn-data')
out_dir = outdir # here('figures', 'venn-plot')
dir.create(out_dir, showWarnings = FALSE)

pattern.str = 'venn.data.*.dat'
print(list.files(data_dir, pattern = pattern.str))

venndata_list <- list()
vennfn_list <- list()
i = 1
for (venfn in list.files(data_dir, pattern = pattern.str)) {
    infn = here('result', 'venn-data', venfn)
    print(infn)
    dt <- read.table(infn)
    base_infn = basename(infn)
    
    if (length(dt$V1) == 7) next
    
    ret = dt$V1
    comb = c("DeepSignal" = ret[1], "Tombo" = ret[2], "Nanopolish" = ret[3],
             "DeepMod" = ret[4], "Megalodon" = ret[5],
             "DeepSignal&Tombo" = ret[6], "DeepSignal&Nanopolish" = ret[7], "DeepSignal&DeepMod" = ret[8],
             "DeepSignal&Megalodon" = ret[9], "Tombo&Nanopolish" = ret[10],
             "Tombo&DeepMod" = ret[11], "Tombo&Megalodon" = ret[12], "Nanopolish&DeepMod" = ret[13],
             "Nanopolish&Megalodon" = ret[14], "DeepMod&Megalodon" = ret[15],
             "DeepSignal&Tombo&Nanopolish" = ret[16], "DeepSignal&Tombo&DeepMod" = ret[17],
             "DeepSignal&Tombo&Megalodon" = ret[18], "DeepSignal&Nanopolish&DeepMod" = ret[19],
             "DeepSignal&Nanopolish&Megalodon" = ret[20], "DeepSignal&DeepMod&Megalodon" = ret[21],
             "Tombo&Nanopolish&DeepMod" = ret[22],
             "Tombo&Nanopolish&Megalodon" = ret[23], "Tombo&DeepMod&Megalodon" = ret[24],
             "Nanopolish&DeepMod&Megalodon" = ret[25],
             "DeepSignal&Tombo&Nanopolish&DeepMod" = ret[26], "DeepSignal&Tombo&Nanopolish&Megalodon" = ret[27],
             "DeepSignal&Tombo&DeepMod&Megalodon" = ret[28], "DeepSignal&Nanopolish&DeepMod&Megalodon" = ret[29],
             "Tombo&Nanopolish&DeepMod&Megalodon" = ret[30],
             "DeepSignal&Tombo&Nanopolish&DeepMod&Megalodon" = ret[31]
    )
    fit1 <- euler(comb, input = 'union', shape = 'circle')
    
    venndata_list[[i]] <- fit1$original.values
    vennfn_list[[i]] <- venfn
    i <- i + 1
}

## TODO: why can not put plot in for loop????


if (dsname == 'NA19240') {
    printi = 4 #NA19240
}

if (dsname == 'APL') {
    printi = 1 #APL
}

if (dsname == 'HL60') {
    printi = 2
}

if (dsname == 'K562') {
    printi = 3
}

outfn = here('figures', 'venn-plot', sprintf('upset.%s.pdf', vennfn_list[[printi]]))
outfn = here('figures', 'venn-plot', sprintf('upset.%s.jpg', vennfn_list[[printi]]))

print(outfn)

upsetTable = fromExpression(venndata_list[[printi]])

# Make marix for ComplexHeatmap usage
# upsetMatrix = make_comb_mat(upsetTable)

dev.new()
dev.list()
# pdf(file = outfn, onefile = FALSE) # or other device

jpg(file = outfn)

upset(upsetTable, scale.intersections = "identity",
      sets.x.label = "Total CpGs", mainbar.y.label = "Intersection CpGs",
      show.numbers = 'no', text.scale = 2.5, number.angles = 0, nintersects = 31, keep.order = TRUE,
      set_size.show = FALSE, set_size.angles = 0, set_size.numbers_size = 6, point.size = 3,
      order.by = "freq")
dev.off()
print(sprintf('save to %s', outfn))

