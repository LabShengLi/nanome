wdir = "/Users/liuya/PycharmProjects/nano-compare/src/rplot_func"
setwd(wdir)

source('utils_plotr.R')

infn = 'performance-results.csv'

# Load data and sort string orders
df <- load.performance.data(infn)

measure.pair.list = list(c('Accuracy', 'ROC_AUC'), c('F1_5mC', 'F1_5C'))

for (measure.pair in measure.pair.list) {
  fig.34a.scatter.plot.performance(df, measure.pair, locations = locations.Singletons)
  fig.34a.scatter.plot.performance(df, measure.pair, locations = locations.Genome)
}


# Plot using functions
perf.measure.list = c('Accuracy', 'AP', 'ROC_AUC', 'F1_5mC', 'Precision_5mC', 'Recall_5mC')

for (perf.measure in perf.measure.list) {
  fig.34a.bar.plot.performance(df, perf.measure, locations = locations.Singletons)
  fig.34a.bar.plot.performance(df, perf.measure, locations = locations.Genome, figsize = c(7, 4))
}



