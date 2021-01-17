library(here)
here()

wdir = here('src', 'rplot_func')  #"/Users/liuya/PycharmProjects/nano-compare/src/rplot_func"

setwd(wdir)

source('utils_plotr.R')

infn = 'performance-results.csv'

# Load data and sort string orders
df <- load.performance.data(infn)

measure.pair.list = list(c('Accuracy', 'Micro.F1'), c('Micro.Precision', 'Micro.Recall'), c('F1_5mC', 'F1_5C'))

for (measure.pair in measure.pair.list) {
  fig.34a.scatter.plot.performance(df, measure.pair, locations = locations.Singletons, scale = 0.65)
  fig.34a.scatter.plot.performance(df, measure.pair, locations = locations.Genome, scale = 0.75)
  #break
}


source('utils_plotr.R')


# Plot using functions
perf.measure.list = c('Accuracy', 'ROC.AUC', 'Micro.F1', 'Macro.F1', 'Average.Precision', 'Recall_5mC')

for (perf.measure in perf.measure.list) {
  fig.34a.bar.plot.performance(df, perf.measure, locations = locations.Singletons, scale = 0.6)
  fig.34a.bar.plot.performance(df, perf.measure, locations = locations.Genome, scale = 0.75)
  #break
}



