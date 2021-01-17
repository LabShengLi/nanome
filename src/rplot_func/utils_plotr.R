library(ggplot2)
library(data.table)

Tool.Order = c('DeepSignal', 'Tombo', 'Nanopolish', 'DeepMod', 'Megalodon')
locations.Singletons = c("Singletons", "Non-singletons", "Discordant", "Concordant")
locations.Genome = c("Genome-wide", "CpG Island", "Promoters", "Exons", "Intergenic", "Introns")
Coord.Order = c(locations.Singletons, locations.Genome)

printf <- function(...) cat(sprintf(...))

load.performance.data <- function(infn) {
  # Load data and sort string orders
  df <- read.csv(file = infn)
  df$Tool <- factor(df$Tool, levels = Tool.Order)
  df$Location <- factor(df$Location, levels = Coord.Order)
  return(df)
}

fig.34a.bar.plot.performance <- function(df, perf.measure = 'Accuracy', locations = locations.Singletons, figsize = c(5.8, 4), scale = 0.65) {
  #Select data in locations
  sel_data = df[df$Location %in% locations,]

  #Plot and save
  outfn = sprintf("figures/fig.34a.bar.%s.%s.png", locations[1], perf.measure)
  p1 <- ggplot(sel_data, aes_string(x = 'Tool', y = perf.measure, fill = 'Tool')) +
    geom_bar(stat = 'identity') +
    facet_grid(Dataset ~ Location) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(strip.text.x = element_text(size = 7))

  ggsave(p1, filename = outfn, scale = scale)
  printf("save to %s\n", outfn)
}

fig.34a.scatter.plot.performance <- function(df, measure.pair, locations, scale = 0.75) {
  # Select data, transfer wide to long format
  sel_data = df[df$Location %in% locations, c('Dataset', 'Tool', 'Location', measure.pair)]
  longdf <- melt(setDT(sel_data), id.vars = c("Dataset", "Tool", "Location"), variable.name = "perf_name")

  # Scatter plot
  p1 <- ggplot(longdf, aes(x = Location, y = value, shape = Tool, color = Location)) +
    geom_point(size = 4) +
    facet_grid(perf_name ~ Dataset, switch = "y") +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5)) +
    ylim(0, 1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10)) +
    guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
    theme(axis.title.y = element_blank()) +
    theme(strip.background.y = element_blank(), strip.placement = "outside")

  outfn = sprintf("figures/fig.34a.scatter.%s.%s.png", locations[1], measure.pair[1])
  ggsave(p1, filename = outfn, scale = scale)
  printf("save to %s\n", outfn)

}