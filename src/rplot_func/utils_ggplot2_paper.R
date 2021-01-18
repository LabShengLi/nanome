library(ggplot2)
library(data.table)

locations.Singletons = c("Singleton", "Nonsingleton", "Discordant", "Concordant")
locations.Genome = c("Genomewide", "CpG Island", "Promoters", "Exons", "Intergenic", "Introns")
Coord.Order = c(locations.Genome, locations.Singletons)
Dataset.Order = c('APL', 'HL60', 'K562')
Tool.Order = c('DeepSignal', 'Tombo', 'Nanopolish', 'DeepMod', 'Megalodon')

measure.pair.list = list(c('Accuracy', 'Micro.F1'), c('Accuracy', 'ROC.AUC'), c('Micro.Precision', 'Micro.Recall'), c('F1_5mC', 'F1_5C'), c('Macro.Precision', 'Macro.Recall'))

measure.list = c('Accuracy', 'ROC.AUC', 'Micro.F1', 'Macro.F1', 'Average.Precision', 'Micro.Precision', 'Micro.Recall', 'Macro.Precision', 'Macro.Recall')

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ToolColorPal <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#F0E442")
ToolShapeList <- c(0, 1, 2, 3, 4, 5)

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
