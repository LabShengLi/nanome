rm(list = ls())
library(here)
library(UpSetR)

data_dir = here('result', 'venn-data')
out_dir = here('figures', 'venn-plot')
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

printi = 2
outfn = here('figures', 'venn-plot', sprintf('upset.%s.pdf', vennfn_list[[printi]]))
print(outfn)
dev.new()
dev.list()
pdf(file = outfn, onefile = FALSE) # or other device

upset(fromExpression(venndata_list[[printi]]), scale.intersections = "identity",
      sets.x.label = "Total CpGs", mainbar.y.label = "Intersection CpGs",
      show.numbers = 'yes', text.scale = 2, number.angles = 30, nintersects = 31, keep.order = TRUE,
      sets.bar.color = "grey", set_size.show = TRUE, set_size.angles = 0, set_size.numbers_size = 6,
      order.by = "freq")
dev.off()
print(sprintf('save to %s', outfn))

