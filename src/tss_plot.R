rm(list = ls())
library(here)
## TSS plot
library(readxl)
library(tidyverse)
library(ggplot2)
library("optparse")

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "/Users/liuya/PycharmProjects/nano-compare/result/NA19240.bin50.outMatrix.for.R.csv.gz",
              help = "Input file"),
  make_option(c("-o", "--output dir"), type = "character", default = getwd(),
              help = "Print little output"),
  make_option(c("--dsname"), type = "character", default = "NA19240",
              help = "Dataset name"),
  make_option(c("--bin-size"), type = "integer", default = 50,
              help = "Bin size used"),
  make_option(c("--bslabel"), type = "character", default = "WGBS",
              help = "WGBS or RRBS")
)

opt = parse_args(OptionParser(option_list = option_list))

print(opt)

flank = 2000
binSize = opt$`bin-size`

numBins = 2 * flank / binSize

width = 8
height = 4

dsname = opt$dsname
infn = opt$input

# infn = '/Users/liuya/PycharmProjects/nano-compare/result/APL.bin50.outMatrix.tsv.gz'
# dsname='APL'

col_types_str = paste0('ciicdc',
                       paste(replicate(numBins * 6, "d"), collapse = ""),
                       collapse = "")

df = read_tsv(infn, col_names = FALSE, col_types = col_types_str, skip = 1)
dim(df)

# Remove first 6 columns, leave only five tools + 1 BS-seq
ndf = df[, 7:(numBins * 6 + 6)]
dim(ndf)

avg_df <- ndf %>%
  summarise_all(mean, na.rm = TRUE)

dim(avg_df)

bsdata = avg_df[, 1:numBins]
deepsignal = avg_df[, (1 + numBins):(numBins * 2)]
tombo = avg_df[, (1 + numBins * 2):(numBins * 3)]
nanopolish = avg_df[, (1 + numBins * 3):(numBins * 4)]
deepmod = avg_df[, (1 + numBins * 4):(numBins * 5)]
megalodon = avg_df[, (1 + numBins * 5):(numBins * 6)]

xvector = seq(-flank, flank - 1, by = binSize)
xlist = xvector
# xlist = rep(xvector, nrow(bsdata))

bslabel = opt$bslabel

list_data <- list("BS-seq" = bsdata, 'DeepSignal' = deepsignal,
                  'Tombo' = tombo, 'Nanopolish' = nanopolish,
                  'DeepMod' = deepmod, 'Megalodon' = megalodon)

list_data

plotdf = tibble()
for (name in names(list_data)) {
  ylist = as.vector(t(list_data[[name]]))
  plotdf1 <- tibble(x = xlist, y = ylist, tool = name)
  plotdf <- bind_rows(plotdf, plotdf1)
}

plotdf = as_tibble(plotdf)
plotdf = transform(plotdf, x = as.numeric(x))
# plotdf = plotdf[1:90000, ]

tools_bsseq = c('DeepSignal', 'Tombo', 'Nanopolish', 'DeepMod', 'Megalodon', bslabel)

plotdf <- plotdf %>%
  mutate(tool = recode(tool, 'BS-seq' = bslabel)) %>%
  mutate(tool = factor(tool, levels = tools_bsseq))

head(plotdf)

dim(plotdf)


print("Start plotting")

color_Pal <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "black", "#0072B2", "#D55E00", "#F0E442")


p1 <- ggplot(plotdf, aes(x = x, y = y, fill = tool, group = tool, color = tool)) +
  # p1 <- ggplot(plotdf, aes(x, y)) +
  # geom_point(size = 0.5) +
  # geom_line() +
  geom_smooth(method = "loess", size = 0.5, se = FALSE, aes(color = tool)) +
  scale_color_manual(values = color_Pal) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # scale_fill_brewer(palette = "Set1") +
  xlim(-2000, 2000) +
  ylim(0, 1) +
  ylab("% methylated") +
  xlab("Binned distance to TSS (bp)")

p1

p2 <- ggplot(plotdf, aes(x = x, y = y, fill = tool, group = tool, color = tool)) +
  geom_point(size = 0.5) +
  geom_line() +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # geom_smooth(method = "loess", size = 1, se = FALSE, aes(color = tool)) +
  scale_color_manual(values = color_Pal) +
  # scale_fill_brewer(palette = "Set1") +
  xlim(-2000, 2000) +
  ylim(0, 1) +
  ylab("% methylated") +
  xlab("Binned distance to TSS (bp)")


print("Start saving")

out_dir = opt$`output dir`
dir.create(out_dir, showWarnings = FALSE)

outfn = sprintf("%s/TSS.%s.binsize%d.smoothed.curves.jpg", out_dir, dsname, binSize)
ggsave(p1, filename = outfn, width = width, height = height, dpi = 600)

outfn2 = sprintf("%s/TSS.%s.binsize%d.lines.curves.jpg", out_dir, dsname, binSize)
ggsave(p2, filename = outfn2, width = width, height = height, dpi = 600)

print(sprintf("save to %s.", outfn))


