library(ggplot2)

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

fig.34a.bar.plot.performance <- function(df, perf.measure = 'Accuracy', locations = locations.Singletons, figsize = c(5.8, 4)) {
  #Select data in locations
  sel_data = df[df$Location %in% locations,]

  #Plot and save
  outfn = sprintf("out/fig.34a.%s.%s.png", locations[1], perf.measure)
  pdf(outfn)
  p1 <- ggplot(sel_data, aes_string(x = 'Tool', y = perf.measure, fill = 'Tool')) +
    geom_bar(stat = 'identity') +
    facet_grid(Dataset ~ Location) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(p1, filename = outfn, width = figsize[1], height = figsize[2])
  printf("save to %s\n", outfn)
}
