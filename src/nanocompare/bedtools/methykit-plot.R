#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#
#BiocManager::install("methylKit")
library(methylKit)

datadir = '/Users/liuya/Downloads/ctimage/12-28'

datadir = '/projects/li-lab/yang/results/2020-12-28/K562_WGBS_Joined'

setwd(datadir)

file.list = as.list(list.files(path = datadir, pattern = 'nocovcutoff', full.names = TRUE))

file.list
# read the files to a methylRawListDB object: myobjDB
# and save in databases in folder methylDB
myobjDB = methRead(file.list,
                   sample.id = list("DeepMod", "DeepSignal", "Nanopolish", "Tombo"),
                   assembly = "hg38",
                   treatment = c(1, 1, 1, 1),
                   context = "CpG",
                   mincov = 4,
                   dbtype = "tabix",
                   dbdir = "methylDB"
)

#getMethylationStats(myobj[[4]], plot = TRUE, both.strands = FALSE)

for (val in seq(1, 4)) {
  outfn = file.path(datadir, sprintf("cov-stats-%d.png", val))

  # png("~/Desktop/Figure4.png", width = 18, height = 21, units = 'cm', res = 300)
  png(outfn, width = 12, height = 14, units = 'cm', res = 400, pointsize = 12)
  getCoverageStats(myobjDB[[val]], plot = TRUE, both.strands = FALSE)
  dev.off()
}