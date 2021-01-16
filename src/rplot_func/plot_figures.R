wdir = "/Users/liuya/PycharmProjects/nano-compare/src/rplot_func"
setwd(wdir)

source('utils_plotr.R')

infn = 'performance-results.csv'

# Load data and sort string orders
df <- load.performance.data(infn)

sel_data = df[df$Location %in% locations.Singletons,]

p1 <- ggplot(sel_data, aes(x = Location, y = Accuracy, shape = Location, color = Tool)) +
  geom_point() +
  facet_grid(~Dataset) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(p1, filename = 'test.png')


# Plot using functions
perf.measure.list = c('Accuracy', 'AP', 'ROC_AUC', 'F1_5mC', 'Precision_5mC', 'Recall_5mC')

for (perf.measure in perf.measure.list) {
  fig.34a.bar.plot.performance(df, perf.measure, locations = locations.Singletons)
  fig.34a.bar.plot.performance(df, perf.measure, locations = locations.Genome, figsize = c(7, 4))
  break
}



