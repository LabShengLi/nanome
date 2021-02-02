library(ggplot2)
library(VennDiagram)
library(data.table)
library(eulerr)
library(tidyverse)
library(here)
library(ggpubr)


locations.Singletons = c("Singleton", "Nonsingleton", "Discordant", "Concordant")
locations.Genome = c("Genomewide", "CpG Island", "Promoters", "Exons", "Intergenic", "Introns")
Coord.Order = c(locations.Genome, locations.Singletons)
Dataset.Order = c('APL', 'HL60', 'K562', 'NA19240')
Tool.Order = c('DeepSignal', 'Tombo', 'Nanopolish', 'DeepMod', 'Megalodon')

measure.pair.list = list(c('Accuracy', 'Micro.F1'), c('Accuracy', 'ROC.AUC'), c('Micro.Precision', 'Micro.Recall'), c('F1_5mC', 'F1_5C'), c('Macro.Precision', 'Macro.Recall'))
measure.list = c('Accuracy', 'ROC.AUC', 'Micro.F1', 'Macro.F1', 'Average.Precision', 'Micro.Precision', 'Micro.Recall', 'Macro.Precision', 'Macro.Recall')
Corr.Perf.List = c('Corr_All', 'Corr_Mix')
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ToolColorPal <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#F0E442")

ToolShapeList <- c(0, 1, 2, 3, 4, 5)

Top3.Tool.Index = c(1, 3, 5)
Top3.Tool.Order = Tool.Order[Top3.Tool.Index]
Top3.ToolColorPal <- ToolColorPal[Top3.Tool.Index]

venn_flist = c('venn.data.HL60.dat', 'venn.data.K562.dat', 'venn.data.APL.dat')


printf <- function(...) cat(sprintf(...))


load.performance.data <- function(infn) {
  # Load data and sort string orders
  df <- read.csv(file = infn)

  df$Location <- as.character(df$Location)
  df[df$Location == 'Genome-wide', 'Location'] <- 'Genomewide'
  df[df$Location == 'Non-singletons', 'Location'] <- 'Nonsingleton'
  df[df$Location == 'Singletons', 'Location'] <- 'Singleton'

  df$Tool <- factor(df$Tool, levels = Tool.Order)
  df$Location <- factor(df$Location, levels = Coord.Order)
  df$Dataset <- factor(df$Dataset, levels = Dataset.Order)
  return(df)
}


fig.34a.bar.plot.performance <- function(df, perf.measure = 'Accuracy', locations = locations.Singletons, bdir, figsize = c(5.8, 4), scale = 0.65) {
  #Select data in locations
  sel_data = df[df$Location %in% locations,]

  #Plot and save
  outfn = sprintf("%s/fig.34a.bar.%s.%s.jpg", bdir, locations[1], perf.measure)
  p1 <- ggplot(sel_data, aes_string(x = 'Tool', y = perf.measure, fill = 'Tool')) +
    geom_bar(stat = 'identity') +
    facet_grid(Dataset ~ Location) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(strip.text.x = element_text(size = 7)) +
    scale_fill_manual(values = ToolColorPal)

  ggsave(p1, filename = outfn, scale = scale)
  printf("save to %s\n", outfn)
}


fig.34a.bar.plot1d.performance <- function(df, perf.measure = 'Accuracy', locations = locations.Singletons, bdir, figsize = c(5.8, 4), scale = 0.65) {
  sel_data = df[df$Location %in% locations, c('Dataset', 'Tool', 'Location', 'Accuracy')]

  ggplot(sel_data, aes_string(x = 'Tool', y = perf.measure, fill = 'Tool')) +
    geom_bar(stat = 'identity') +
    facet_grid(~Dataset + Location) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(strip.text.x = element_text(size = 7)) +
    scale_fill_manual(values = cbPalette)

  outfn = sprintf("%s/bar1d.%s.%s.jpg", bdir, perf.measure, locations[1])
  ggsave(filename = outfn, width = 10, height = 4, dpi = 600, limitsize = FALSE)

}


fig.34a.scatter.plot.performance <- function(df, measure.pair, locations, bdir, scale = 0.75) {
  # Select data, transfer wide to long format
  sel_data = df[df$Location %in% locations, c('Dataset', 'Tool', 'Location', measure.pair)]
  longdf <- melt(setDT(sel_data), id.vars = c("Dataset", "Tool", "Location"), variable.name = "perf_name")

  # Scatter plot
  p1 <- ggplot(longdf, aes(x = Location, y = value, shape = Tool, color = Location)) +
    geom_point(size = 4) +
    facet_grid(perf_name ~ Dataset, switch = "y") +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5)) +
    #scale_fill_manual(values = cbPalette) +
    scale_color_manual(values = cbPalette) +
    ylim(0, 1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10)) +
    guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
    theme(axis.title.y = element_blank()) +
    theme(strip.background.y = element_blank(), strip.placement = "outside")

  outfn = sprintf("%s/fig.34a.scatter.%s.%s.jpg", bdir, locations[1], measure.pair[1])
  ggsave(p1, filename = outfn, scale = scale)
  printf("save to %s\n", outfn)

}


fig.34a.line.plot.performance <- function(df, perf.measure, locations, bdir, scale = 1) {
  sel_df = df[df$Location %in% locations, c('Dataset', 'Tool', 'Location', perf.measure)]

  ggplot(data = sel_df, mapping = aes_string(x = 'Location', y = perf.measure, group = 'Tool')) +
    facet_grid(~Dataset) +
    geom_point(aes(shape = Tool, color = Tool), size = 5) +
    geom_line(aes(linetype = Tool, color = Tool)) +
    ylim(0, 1) +
    scale_shape_manual(values = ToolShapeList) +
    scale_color_manual(values = ToolColorPal) +
    scale_linetype_manual(values = c("dashed", "dotted", "twodash", "dotdash", "longdash", "twodash")) +
    theme(legend.position = "top") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)) +
    theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10)) +
    theme(strip.text.x = element_text(size = 12))

  outfn = sprintf("%s/fig.34a.line.%s.%s.jpg", bdir, perf.measure, locations[1])
  ggsave(filename = outfn, width = 6.5, height = 4, scale = scale, dpi = 600)
  printf("save to %s\n", outfn)
}


fig.34a.bar.dataset.location.performance <- function(df, dsname, perf.pair, bdir, scale = 1) {
  sel_data = df[df$Dataset == dsname, c('Dataset', 'Tool', 'Location', perf.pair)]
  longdf <- melt(setDT(sel_data), id.vars = c("Dataset", "Tool", "Location"), variable.name = "perf_name")

  ggplot(longdf, aes_string(x = 'Tool', y = 'value', fill = 'Tool')) +
    geom_bar(stat = 'identity') +
    facet_grid(perf_name ~ Location, switch = "y") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)) +
    theme(strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 12)) +
    scale_fill_manual(values = ToolColorPal) +
    theme(axis.title.y = element_blank()) +
    theme(strip.background.y = element_blank(), strip.placement = "outside") +
    theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10))

  outfn = sprintf("%s/fig.34a.bar.dataset.%s.%s.jpg", bdir, dsname, perf.pair[1])
  ggsave(filename = outfn, width = 10, height = 4, scale = scale, dpi = 600)
  printf("save to %s\n", outfn)
}


fig.6a.venn.plot.set5 <- function(ret, outfn) {
  graphics.off()
  venn.plot <- draw.quintuple.venn(
    area1 = ret[1],
    area2 = ret[2],
    area3 = ret[3],
    area4 = ret[4],
    area5 = ret[5],
    n12 = ret[6],
    n13 = ret[7],
    n14 = ret[8],
    n15 = ret[9],
    n23 = ret[10],
    n24 = ret[11],
    n25 = ret[12],
    n34 = ret[13],
    n35 = ret[14],
    n45 = ret[15],
    n123 = ret[16],
    n124 = ret[17],
    n125 = ret[18],
    n134 = ret[19],
    n135 = ret[20],
    n145 = ret[21],
    n234 = ret[22],
    n235 = ret[23],
    n245 = ret[24],
    n345 = ret[25],
    n1234 = ret[26],
    n1235 = ret[27],
    n1245 = ret[28],
    n1345 = ret[29],
    n2345 = ret[30],
    n12345 = ret[31],
    category = Tool.Order,
    fill = ToolColorPal[1:5],
    cat.col = ToolColorPal[1:5],
    cat.pos = c(0, 2, 8, 3, 11),
    cat.dist = c(0.20, 0.22, -0.17, -0.18, 0.22),
    cat.cex = 1.8,
    margin = 0.05,
    cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
            1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
    ind = TRUE
  );

  ggsave(venn.plot, file = outfn, width = 5, height = 5, dpi = 600)
  printf('save to %s', outfn)
}


fig.6a.venn.plot.set3 <- function(ret, outfn) {
  graphics.off()
  venn.plot <- draw.triple.venn(
    area1 = ret[1],
    area2 = ret[2],
    area3 = ret[3],
    n12 = ret[4],
    n23 = ret[6],
    n13 = ret[5],
    n123 = ret[7],
    category = Top3.Tool.Order,
    fill = Top3.ToolColorPal,
    lty = "blank",
    cex = 2,
    cat.cex = c(3, 3, 3),
    cat.col = Top3.ToolColorPal,
    ind = FALSE,
    cat.pos = c(-20, 20, 180),
    print.mode = c("raw", "percent"),
    fontfamily = rep("Times New Roman", 7),
    cat.fontfamily = rep("Times New Roman", 3),
    euler.d = TRUE, scaled = TRUE
  );

  #print("Add title")
  #require(gridExtra)
  #grid.arrange(gTree(children = venn.plot), top = basename(outfn))

  #print(venn.plot)

  ggsave(venn.plot, file = outfn, dpi = 600)
  printf('save to %s', outfn)
}


fig.6b.euller.plot.set3 <- function(ret, outfn) {
  graphics.off()

  fit1 <- euler(c("A" = ret[1], "B" = ret[2], "C" = ret[3],
                  "A&B" = ret[4], "A&C" = ret[5], "B&C" = ret[6],
                  "A&B&C" = ret[7]), input = 'union', shape = 'circle')
  gp <- plot(fit1,
             quantities = list(type = c("counts"),
                               font = 1, round = 2, cex = 0.6),
             labels = identical(legend, FALSE),
             fills = Top3.ToolColorPal,
             alpha = 0.5,
             fill_opacity = 0.5,
             edges = list(lty = 1),
             #main = basename(outfn),
             counts = TRUE,
             legend = list(labels = Top3.Tool.Order, cex = 0.5,
                           colors = Top3.ToolColorPal, alpha = 1, side = 'bottom', ncol = 3, nrow = 1, hgap = 0),
             adjust_labels = TRUE,
             pal = Top3.ToolColorPal,
  )
  ggsave(gp, file = outfn, dpi = 600, width = 3, height = 2)
}


fig.5d.violin.corr.performance <- function(df, corr.perf, bdir, scale = 1) {
  graphics.off()
  ggplot(df, aes_string(x = 'Tool', y = corr.perf, fill = 'Tool')) +
    geom_violin() +
    scale_fill_manual(values = ToolColorPal) +
    geom_jitter(shape = 16, position = position_jitter(0.2), size = 0.5) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1.2, size = 10)) +
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8)) +
    theme(legend.position = "top") +
    guides(fill = guide_legend(ncol = 3, nrow = 2, byrow = TRUE))

  #scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  #NULL

  outfn = sprintf("%s/fig.34a.violin.corr.%s.jpg", bdir, corr.perf)
  ggsave(filename = outfn, width = 3.5, height = 4, scale = scale, dpi = 600)
  printf("save to %s\n", outfn)
}


fig.34b.box.location.performance <- function(df, perf.measure, locations, bdir, scale = 1) {
  sel_data = df[df$Location %in% locations, c('Dataset', 'Tool', 'Location', perf.measure)]

  ggplot(sel_data, aes_string(x = 'Tool', y = perf.measure, fill = 'Tool')) +
    #geom_violin() +
    geom_boxplot() +
    facet_grid(~Location) +
    scale_fill_manual(values = ToolColorPal) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10))

  outfn = sprintf("%s/fig.34a.box.perfmeasure.%s.%s.jpg", bdir, perf.measure, locations[1])
  ggsave(filename = outfn, width = 7.5, height = 3, scale = scale, dpi = 600)
  printf("save to %s\n", outfn)

}


fig.5c.running.resource.bar.plot <- function(bdir) {
  infn1 = here('result', 'recalculate.running.summary.csv')
  run.table1 <- read_csv(infn1)

  infn2 = here('result', 'recalculate.running.summary.Megalodon.csv')
  run.table2 <- read_csv(infn2)

  infn3 = here('result', 'recalculate.running.summary.na19240.csv')
  run.table3 <- read_csv(infn3)

  run.table = bind_rows(run.table1, run.table2, run.table3)


  run.table$rt <- run.table$running.time.seconds / run.table$fast5
  run.table$mem <- run.table$mem.usage.gb / run.table$fast5 * 1000
  run.table$tool <- factor(run.table$tool, levels = Tool.Order)
  run.table$dsname <- factor(run.table$dsname, levels = Dataset.Order)

  run.mean.table <- run.table %>%
    group_by(dsname, tool) %>%
    summarise(mean.rt = mean(rt), mean.mem = mean(mem))

  outfn = here('result', 'mean.runnning.summary.csv')
  write_csv(run.mean.table, outfn)

  g1 <- run.table %>%
    drop_na(tool) %>%
    ggplot(mapping = aes(x = dsname, y = rt, fill = tool)) +
    geom_bar(stat = "summary", fun = mean, width = 0.8, position = position_dodge()) +
    scale_fill_manual(values = ToolColorPal) +
    theme_classic() +
    ylab("Running time(seconds) per Fast5 file") +
    xlab("Dataset") +
    coord_flip() +
    labs(fill = 'Tool') +
    theme(text = element_text(size = 12))

  #g1

  g2 <- run.table %>%
    drop_na(tool) %>%
    ggplot(mapping = aes(x = dsname, y = mem, fill = tool)) +
    geom_bar(stat = "summary", fun = mean, width = 0.8, position = position_dodge()) +
    #scale_x_discrete(limits = levels(Tool.Order)) +
    coord_flip() +
    scale_fill_manual(values = ToolColorPal) +
    theme_classic() +
    ylab("Memory usage(MB) per Fast5 file") +
    xlab("Dataset") +
    labs(fill = 'Tool') +
    theme(text = element_text(size = 12))

  #g2

  gg <- ggarrange(g1, g2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")

  #gg

  outfn = sprintf("%s/running.resource.bar.plot.jpg", bdir)
  ggsave(gg, filename = outfn, width = 8, height = 3, dpi = 600)
  printf("save to %s\n", outfn)
}