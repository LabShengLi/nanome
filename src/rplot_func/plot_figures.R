wdir = "/Users/liuya/PycharmProjects/nano-compare/src/rplot_func"
setwd(wdir)

source('utils_plotr.R')

infn = 'performance-results.csv'

# Load data and sort string orders
df <- load.performance.data(infn)

# Plot using functions
perf.measure.list = c('Accuracy', 'AP', 'ROC_AUC', 'F1_5mC', 'Precision_5mC', 'Recall_5mC')

for (perf.measure in perf.measure.list) {
  fig.34a.bar.plot.performance(df, perf.measure, locations = locations.Singletons)
  fig.34a.bar.plot.performance(df, perf.measure, locations = locations.Genome)
}



