library(ggplot2)
library(VennDiagram)
library(data.table)
library(eulerr)
library(tidyverse)
library(here)
library(ggpubr)

infn.perf = here('result', 'performance-results.csv')
infn.corr = here('result', 'All.corrdata.coe.pvalue.each.regions.xlsx')

locations.Singletons = c("Genome-wide", "Singletons", "Non-singletons", "Concordant", "Discordant")
locations.Regions = c("Promoters", "Exons", "Introns", "Intergenic", "CpG island", "CpG shore", "CpG shelf")
Coord.Order = c(locations.Regions, locations.Singletons)
Dataset.Order = c('HL60', 'K562', 'APL', 'NA19240')
Tool.Order = c('DeepSignal', 'Tombo', 'Nanopolish', 'DeepMod', 'Megalodon', 'Joined', 'Union')
Tool.Order.show.SingleRead = c('DeepSignal', 'Tombo', 'Nanopolish*', 'DeepMod', 'Megalodon*', 'Joined', 'Union')


measure.pair.list = list(c('Accuracy', 'Micro.F1'), c('Accuracy', 'ROC.AUC'), c('Micro.Precision', 'Micro.Recall'), c('F1_5mC', 'F1_5C'), c('Macro.Precision', 'Macro.Recall'))
#measure.list = c('Accuracy', 'ROC.AUC', 'Micro.F1', 'Macro.F1', 'Average.Precision', 'Micro.Precision', 'Micro.Recall', 'Macro.Precision', 'Macro.Recall')
#measure.list = c('Accuracy', 'ROC.AUC', 'Macro.F1', 'Average.Precision', 'Macro.Precision', 'Macro.Recall')
measure.list = c('Accuracy', 'ROC.AUC', 'Macro.F1')

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


export.table.s3.xlsx <- function() {
  locations.Singletons = c("Genome-wide", "Singletons", "Non-singletons", "Concordant", "Discordant")
  locations.Genome = c("Promoters", "Exons", "Introns", "Intergenic", "CpG island", "CpG shore", "CpG shelf")

  library(readxl)
  library(tidyverse)
  out_dir = here('figures')
  df <- read_csv(infn.perf)

  selected.columns = c('referenceCpGs', 'mCsites', 'Csites')
  outdf <- df %>%
    mutate(Location = recode(Location, 'CpG Island' = 'CpG island', 'CpG Shores' = 'CpG shore', 'CpG Shelves' = 'CpG shelf')) %>%
    mutate(Dataset = factor(Dataset, levels = Dataset.Order),
           Location = factor(Location, levels = c(locations.Singletons, locations.Genome))) %>%
    arrange(desc(Dataset), Location, desc(Accuracy)) %>%
    select(seq(2, 10), selected.columns) %>%
    drop_na()

  outfn = here("figures", "Table.S3.csv")
  write_csv(outdf, outfn)
  print(sprintf("save to %s", outfn))

}


export.table.s4.xlsx <- function() {
  library(readxl)
  library(tidyverse)

  df = read_excel(infn.corr)

  locations.Singletons = c("Genome-wide", "Singletons", "Non-singletons", "Concordant", "Discordant")
  locations.Genome = c("Promoters", "Exons", "Introns", "Intergenic", "CpG island", "CpG shore", "CpG shelf")

  outdf <- df %>%
    mutate(Location = recode(Location, 'CpG Island' = 'CpG island', 'CpG Shores' = 'CpG shore', 'CpG Shelves' = 'CpG shelf')) %>%
    mutate(dsname = factor(dsname, levels = Dataset.Order),
           Location = factor(Location, levels = c(locations.Singletons, locations.Genome))) %>%
    arrange(desc(dsname), Location, desc(COE)) %>%
    drop_na() %>%
    select(2:ncol(df))

  outfn = here("figures", "Table.S4.csv")
  write_csv(outdf, outfn)


}


load.performance.data <- function() {
  library(readxl)
  library(tidyverse)

  # Load data and sort string orders
  df <- read.csv(infn.perf)

  outdf <- df %>%
    mutate(Location = recode(Location, 'CpG Island' = 'CpG island', 'CpG Shores' = 'CpG shore', 'CpG Shelves' = 'CpG shelf')) %>%
    mutate(Dataset = factor(Dataset, levels = Dataset.Order),
           Location = factor(Location, levels = c(locations.Singletons, locations.Regions)),
           Tool = factor(Tool, levels = Tool.Order)
    ) %>%
    drop_na()

  return(outdf)
}


fig.s34a.line.plot.performance <- function(perf.measure, locations, bdir, scale = 1) {
  df <- load.performance.data()

  sel_df = df[df$Location %in% locations, c('Dataset', 'Tool', 'Location', perf.measure)]

  out.y.label = perf.measure
  if (perf.measure == 'ROC.AUC') {
    out.y.label = 'ROC AUC'
  }

  if (perf.measure == 'Average.Precision') {
    out.y.label = 'Average Precision'
  }

  if (perf.measure == 'Macro.F1') {
    out.y.label = 'F1'
  }
  p <- ggplot(data = sel_df, mapping = aes_string(x = 'Location', y = perf.measure, group = 'Tool')) +
    facet_grid(~Dataset) +
    geom_point(aes(shape = Tool, color = Tool), size = 5) +
    geom_line(aes(linetype = Tool, color = Tool)) +
    ylim(0, 1) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      #panel.grid.major = element_blank(),
      #panel.border = element_blank(),
      #axis.ticks = element_line(size = 0),
      #panel.grid.minor.y = element_blank(),
      #panel.grid.major.y = element_blank(),
      strip.text.x = element_text(size = 8, face = "bold"),
      #axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1, size=0.8)
      #legend.position = "none"
    ) +
    scale_shape_manual(values = ToolShapeList) +
    scale_color_manual(values = ToolColorPal) +
    scale_linetype_manual(values = c("dashed", "dotted", "twodash", "dotdash", "longdash", "twodash")) +
    theme(legend.position = "top") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)) +
    theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10)) +
    theme(strip.text.x = element_text(size = 12)) +
    labs(y = out.y.label)

  outfn = sprintf("%s/fig.34a.line.%s.%s.jpg", bdir, perf.measure, locations[2])
  ggsave(p, filename = outfn, width = 6.5, height = 4, scale = scale, dpi = 600)
  printf("save to %s\n", outfn)
}


fig.34a.box.location.performance <- function(perf.measure, locations, bdir, scale = 1) {
  df <- load.performance.data()

  sel_data = df[df$Location %in% locations, c('Dataset', 'Tool', 'Location', perf.measure)]
  out.y.label = perf.measure
  if (perf.measure == 'ROC.AUC') {
    out.y.label = 'ROC AUC'
  }

  if (perf.measure == 'Average.Precision') {
    out.y.label = 'Average Precision'
  }

  if (perf.measure == 'Macro.F1') {
    out.y.label = 'F1'
  }

  p <- ggplot(sel_data, aes_string(x = 'Tool', y = perf.measure, fill = 'Tool')) +
    #geom_violin() +
    geom_boxplot(outlier.size = 0.2) +
    facet_grid(~Location) +
    scale_fill_manual(values = ToolColorPal) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      #panel.grid.major = element_blank(),
      #panel.border = element_blank(),
      #axis.ticks = element_line(size = 0),
      #panel.grid.minor.y = element_blank(),
      #panel.grid.major.y = element_blank(),
      strip.text.x = element_text(size = 8.5, face = "bold"),
      #strip.text.x = element_text(size = 10),

      panel.border = element_rect(colour = "black"),
      #axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1, size=0.8)
      #legend.position = "none"
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
    ylim(0, 1) +
    labs(y = out.y.label)

  outfn = sprintf("%s/fig.34a.box.perfmeasure.%s.%s.jpg", bdir, perf.measure, locations[2])
  ggsave(p, filename = outfn, width = 7.5, height = 2.5, dpi = 600)
  printf("save to %s\n", outfn)

}


fig.5cd.venn.plot.set5 <- function(ret, outfn) {
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
    category = Tool.Order[1:5],
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


fig.5cd.euller.plot.set3 <- function(ret, outfn) {
  graphics.off()

  comb = c("A" = ret[1], "B" = ret[2], "C" = ret[3],
           "A&B" = ret[4], "A&C" = ret[5], "B&C" = ret[6],
           "A&B&C" = ret[7])

  print(comb)

  fit1 <- euler(comb, input = 'union', shape = 'circle')
  #fit1 <- euler(comb, input = 'union', shape = 'ellipse')

  #fit1 <- venn(comb, input = 'union')

  print(fit1)

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


fig.6a.bar.plot.tools.sites.all.datasets <- function() {
  data_dir = here('result')
  out_dir = here('figures', 'bar-plot')
  dir.create(out_dir, showWarnings = FALSE)
  pattern.str = '*-summary-bgtruth-tools-bsCov5-minCov3.csv'

  totaldt = tibble()
  for (fn in list.files(data_dir, pattern = pattern.str)) {
    infn = here('result', fn)
    basename_infn = basename(infn)

    pos = str_locate(basename_infn, "_")[1]
    dsname = substr(basename_infn, 1, pos - 1)
    dt <- read_csv(infn)
    seldt <- dt[c(1, 2, 3, 4, 5, 6, 7), c(1, 5)]
    names(seldt) <- c('Tool', 'Sites')
    seldt$Tool <- factor(seldt$Tool, levels = Tool.Order)
    seldt$dsname = dsname

    totaldt <- bind_rows(totaldt, seldt)

  }

  plot_totaldt <- totaldt %>%
    mutate(dsname = factor(dsname, levels = Dataset.Order[3:4])) %>%
    mutate(Tool = factor(Tool, levels = Tool.Order)) %>%
    drop_na()

  # totaldt$dsname <- factor(totaldt$dsname, levels = Dataset.Order)
  # totaldt$Tool <- factor(totaldt$Tool, levels = Tool.Order)

  library(scales)

  p1 <- ggplot(data = plot_totaldt, aes(x = Tool, y = Sites, fill = Tool)) +
    geom_bar(stat = 'identity') +
    facet_wrap(~dsname, ncol = 4) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      # panel.grid.major = element_line(colour = "grey80"),
      #panel.border = element_blank(),
      # panel.border = element_rect(colour = "black"),
      #axis.ticks = element_line(size = 0),
      #panel.grid.minor.y = element_blank(),
      #panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.text.x = element_text(size = 12, face = "bold"),
      #axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1, size=0.8)
      #legend.position = "none"
    ) +
    theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1)) +
    theme(strip.text.x = element_text(size = 12)) +
    scale_fill_manual(values = ToolColorPal) +
    ylab("Number of CpGs(cov>=3)") +
    ylim(0, 58000000) +
    scale_y_continuous(labels = comma)


  #+
  #theme(legend.title = element_blank()) +
  # scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  #               labels = scales::trans_format("log10", scales::math_format(10^.x)))
  p1

  outfn = sprintf("%s/fig.5cd.bar.sites.of.tools.all.data.facet.plot.jpg", out_dir)
  ggsave(p1, filename = outfn, width = 5, height = 3, dpi = 600)
  printf("save to %s\n", outfn)
}


fig.6.running.resource.bar.plot <- function() {
  outdir = here('figures')
  infn1 = here('result', 'recalculate.running.summary.csv')
  run.table1 <- read_csv(infn1)

  infn2 = here('result', 'recalculate.running.summary.Megalodon.csv')
  run.table2 <- read_csv(infn2)

  infn3 = here('result', 'recalculate.running.summary.na19240.csv')
  run.table3 <- read_csv(infn3)

  run.table = bind_rows(run.table1, run.table2, run.table3)

  run.table <- run.table %>%
    group_by(dsname, tool) %>%
    summarise(rt = sum(running.time.seconds), mem = max(mem.usage.gb)) %>%
    mutate(rt = rt / 60 / 60)

  run.table$tool <- factor(run.table$tool, levels = Tool.Order)
  run.table$dsname <- factor(run.table$dsname, levels = Dataset.Order)

  run.table <- run.table %>%
    drop_na(tool) %>%
    arrange(dsname, tool)

  outfn = here('result', 'total.runnning.summary.csv')
  write_csv(run.table, outfn)

  g1 <- run.table %>%
    drop_na(tool) %>%
    ggplot(mapping = aes(x = dsname, y = rt, fill = tool)) +
    geom_bar(stat = "summary", fun = mean, width = 0.8, position = position_dodge()) +
    scale_fill_manual(values = ToolColorPal) +
    theme_classic() +
    ylab("Running time (hours)") +
    scale_y_log10() +
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
    ylab("Memory usage (GB)") +
    xlab("") +
    labs(fill = 'Tool') +
    theme(text = element_text(size = 12))

  #g2

  gg <- ggarrange(g1, g2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")

  gg

  outfn = sprintf("%s/fig.6.running.resource.bar.plot.jpg", outdir)
  ggsave(gg, filename = outfn, width = 6, height = 3, dpi = 600)
  printf("save to %s\n", outfn)
}


fig.7b1.running.resource.bar.plot <- function() {
  library(readxl)
  library(tidyverse)

  outdir = here('figures')
  infn = here('result', 'running.logs.total.resource.usage.five.tools.xlsx')
  run.table <- read_excel(infn)

  outdf <- run.table %>%
    mutate(dsname = factor(dsname, levels = Dataset.Order),
           tool = factor(tool, levels = Tool.Order[1:5])
    ) %>%
    arrange(dsname, `Job Wall-clock Time`)

  outfn = here('figures', 'table.s6.total.runnning.summary.new.sorted.csv')
  write_csv(outdf, outfn)

  colnames(outdf) = c('dsname', 'tool', 'cpu.time', 'wall.time', 'mem.usage')

  g1 <- outdf %>%
    drop_na(tool) %>%
    ggplot(mapping = aes(x = cpu.time, y = dsname, fill = tool)) +
    geom_bar(stat = "identity", width = 0.8, position = position_dodge()) +
    scale_fill_manual(values = ToolColorPal) +
    theme_classic() +
    xlab("CPU Utilized Time (hours)") +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides = 'b') +
    ylab("Dataset") +
    labs(fill = 'Tool') +
    theme(text = element_text(size = 12))


  g2 <- outdf %>%
    drop_na(tool) %>%
    ggplot(mapping = aes(x = dsname, y = mem.usage, fill = tool)) +
    geom_bar(stat = "summary", fun = mean, width = 0.8, position = position_dodge()) +
    #scale_x_discrete(limits = levels(Tool.Order)) +
    coord_flip() +
    scale_fill_manual(values = ToolColorPal) +
    theme_classic() +
    ylab("Peak memory (GB)") +
    xlab("") +
    labs(fill = 'Tool') +
    theme(text = element_text(size = 12))

  #g2

  gg <- ggarrange(g1, g2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")

  gg

  outfn = sprintf("%s/fig.7.running.resource.bar.plot.jpg", outdir)
  ggsave(gg, filename = outfn, width = 6.5, height = 3.5, dpi = 600)
  printf("save to %s\n", outfn)

}


fig.7b2.running.resource.bar.plot <- function() {
  library(readxl)
  library(tidyverse)

  outdir = here('figures')
  infn = here('result', 'benchmarking.log.formated.table.step2.all.tools.csv')
  df <- read_csv(infn)

  outdf <- df %>%
    mutate(tool = factor(tool, levels = Tool.Order[1:5]))

  p1 <- ggplot(data = outdf, aes(x = reads, y = realtime, group = tool, colour = tool)) +
    geom_line() +
    geom_point(aes(shape = tool), size = 5) +
    # geom_smooth(method = "loess") +
    scale_color_manual(values = ToolColorPal) +
    scale_shape_manual(values = ToolShapeList) +
    theme_classic() +
    theme(text = element_text(size = 12)) +
    ylab("Wall-clock time (seconds)") +
    xlab("Number of reads") #+
  # scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  #              labels = scales::trans_format("log10", scales::math_format(10^.x)))


  p1


  p2 <- ggplot(data = outdf, aes(x = reads, y = peak_rss, group = tool, colour = tool)) +
    geom_line() +
    geom_point(aes(shape = tool), size = 5) +
    # geom_smooth(method = "loess") +
    scale_color_manual(values = ToolColorPal) +
    scale_shape_manual(values = ToolShapeList) +
    theme_classic() +
    theme(text = element_text(size = 12)) +
    ylab("Peak memory (GB)") +
    xlab("Number of reads") +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides = 'l')


  p2


  outfn = sprintf("%s/fig.7B1.benchmarking.running.time.jpg", outdir)
  ggsave(p1, filename = outfn, width = 4, height = 3, dpi = 600)
  printf("save to %s\n", outfn)


  outfn = sprintf("%s/fig.7B2.benchmarking.peak.memory.jpg", outdir)
  ggsave(p2, filename = outfn, width = 4, height = 3, dpi = 600)
  printf("save to %s\n", outfn)

}


fig.5b.sorted.bar.plot.coe.in.each.region <- function() {
  library(readxl)
  library(tidyverse)
  out_dir = here('figures', 'bar-plot')
  #infn = here('result', 'All.corrdata.coe.pvalue.each.regions.xlsx')
  df = read_excel(infn.corr)

  outdf <- df %>%
    mutate(Location = recode(Location, 'CpG Island' = 'CpG island', 'CpG Shores' = 'CpG shore', 'CpG Shelves' = 'CpG shelf')) %>%
    mutate(dsname = factor(dsname, levels = Dataset.Order),
           Tool = factor(Tool, levels = Tool.Order[1:5]),
           Location = factor(Location, levels = c(locations.Singletons[2:5], locations.Singletons[1:1], locations.Regions[1:7])),

    ) #%>%

  ## Bar plot using facet grid, and order Tools by values
  p1 <- outdf %>%
    drop_na(Tool, Location) %>%
    filter(dsname == 'NA19240') %>%
    mutate(ReorderTool = tidytext::reorder_within(Tool, COE, Location)) %>%
    ggplot(aes(ReorderTool, COE, fill = Tool)) +
    geom_col() +
    coord_flip() +
    tidytext::scale_x_reordered() +
    facet_wrap(vars(Location), scales = "free_y", ncol = 3, dir = 'v') +
    scale_fill_manual(values = ToolColorPal[1:5]) +
    labs(y = "Pearson correlation coefficient", x = NULL) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid.major = element_line(colour = "grey80"),
      panel.border = element_blank(),
      axis.ticks = element_line(size = 0),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      strip.text.x = element_text(size = 12, face = "bold"),
      #axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1, size=0.8)
      #legend.position = "none"
    ) +
    ylim(min(df$COE), 1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
    theme(axis.text.y = element_text(size = 10)) +
    theme(legend.position = "right", legend.text = element_text(size = 10))

  outfn = sprintf("%s/fig.5b.bar.plot.coe.in.each.regions.jpg", out_dir)
  ggsave(p1, filename = outfn, width = 9, height = 7.5, dpi = 600)
  printf("save to %s\n", outfn)
}


fig.s1.pie.plot.singletons.nonsingletons.raw.fast5 <- function() {
  library(readxl)
  library(tidyverse)
  infn = here('result', 'raw.fast5.reads.cpg.coverage.across.regions.cutoff3.xlsx')
  df = read_excel(infn)

  meltdf = df[, c(2, 6, 7)]

  colnames(meltdf) <- c('dsname', 'Singletons', 'Nonsingletons')

  meltdf$dsname <- factor(meltdf$dsname, levels = Dataset.Order)

  row_sum = rowSums(select(meltdf, -dsname))

  #meltdf <- meltdf %>% select(-dsname) %>% mutate_all(~ ./row_sum) %>% head()

  meltdf <- meltdf %>%
    mutate(row_sum = rowSums(select(., 2:3))) %>%
    mutate_at(2:3, ~. / row_sum) %>%
    select(-row_sum) %>%
    head()

  meltdf <- meltdf %>%
    pivot_longer(!dsname, names_to = "type", values_to = "count") %>%
    mutate(label = paste0(round(count * 100, 1), "%"))

  meltdf$type <- factor(meltdf$type, levels = c('Singletons', 'Nonsingletons'))


  p1 <- ggplot(data = meltdf, aes(x = "", y = count, fill = type)) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      axis.ticks = element_line(size = 0),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      strip.text.x = element_text(size = 12, face = "bold"),
      #axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1, size=0.8)
      #legend.position = "none"
    ) +
    geom_bar(stat = "identity", position = position_fill()) +
    geom_text(aes(label = label), position = position_fill(vjust = 0.5)) +
    coord_polar(theta = "y") +
    facet_wrap(~dsname) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    theme(legend.position = 'bottom', legend.title = element_blank()) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    scale_fill_brewer(palette = "Dark2")

  p1
  out_dir = here('figures')
  outfn = sprintf("%s/pie.chart.singleton.nonsingleton.raw.fast5.read.cov3.jpg", out_dir)
  ggsave(p1, filename = outfn, width = 3, height = 4, dpi = 600)
  printf("save to %s\n", outfn)
}


pie.plot.supf1 <- function() {
  # Prepare a color palette. Here with R color brewer:
  library(RColorBrewer)
  myPalette <- brewer.pal(5, "Set1")

  piedatafn = here('result', 'dataset.singleton.vs.non-singleton.csv')
  piedata <- read.csv(file = piedatafn)

  out_dir = here('figures', 'pie-plot')
  dir.create(out_dir, showWarnings = FALSE)

  for (i in 1:nrow(piedata)) {
    row <- piedata[i,]
    dsname = row$Dataset
    total.singletons = row$Singleton.5mc + row$Singleton.5C
    total.concordant = row$Concordant.5mC + row$Concordant.5C
    total.discordant = row$Discordant.5mC + row$Discordant.5C

    data <- c(total.singletons, total.concordant + total.discordant)
    pct <- round(data / sum(data) * 100)
    lbls <- c('Singletons', 'Non-Singletons')
    #lbls <- paste(lbls, ":\n", sep = "") # ad % to labels
    lbls <- paste(lbls, pct) # add percents to labels
    lbls <- paste(lbls, "%", sep = "") # ad % to labels

    outfn = sprintf("%s/pie.plot.%s.jpg", out_dir, dsname)
    jpeg(filename = outfn, width = 3, height = 3, units = "in", res = 600)

    # You can change the border of each area with the classical parameters:
    pie(data, labels = lbls,
        border = "white", col = myPalette, main = dsname, cex = 0.5)

    dev.off()

    #ggsave(filename = outfn, dpi = 600)
    printf("save to %s\n", outfn)

    #break
  }

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


fig.5d.violin.corr.performance <- function(df, bdir) {
  tdf <- as_tibble(df)
  seldf <- tdf %>% select(Corr.Perf.List, 'Tool')

  longdf <- seldf %>%
    pivot_longer(!Tool, names_to = "value_name", values_to = "Pearson_COE")

  longdf <- longdf[longdf$Pearson_COE >= 0,]

  longdf <- longdf %>% drop_na(Pearson_COE)

  graphics.off()
  ggplot(longdf, aes_string(x = 'Tool', y = 'Pearson_COE', fill = 'Tool')) +
    geom_violin() +
    facet_grid(. ~ value_name) +
    scale_fill_manual(values = ToolColorPal) +
    geom_jitter(shape = 16, position = position_jitter(0.2), size = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1.2, size = 10)) +
    theme(legend.position = "top") +
    theme(legend.title = element_blank())

  outfn = sprintf("%s/fig.5d.violin.corr.violin.jpg", bdir)
  ggsave(filename = outfn, width = 5.5, height = 3, dpi = 600)
  printf("save to %s\n", outfn)
}
