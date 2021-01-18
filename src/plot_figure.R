rm(list = ls())

library(here)
here()

getwd()
#wdir = here('src')  #"/Users/liuya/PycharmProjects/nano-compare/src/rplot_func"
#setwd(wdir)

source(here('src', 'rplot_func', 'utils_ggplot2_paper.R'))

outdir = here('figures')

# Load data and sort string orders
infn = here('result', 'performance-results.csv')
df <- load.performance.data(infn)


#### Test


sel_data = df[df$Location %in% locations.Singletons, c('Dataset', 'Tool', 'Location', 'Accuracy')]

ggplot(sel_data, aes_string(x = 'Tool', y = 'Accuracy', fill = 'Tool')) +
    geom_bar(stat = 'identity') +
    facet_grid( ~ Dataset+Location) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(strip.text.x = element_text(size = 7))

outfn = sprintf("%s/test.pdf", outdir)
ggsave(filename = outfn, width=9, height=4, limitsize =FALSE)


# Scatter plot
measure.pair.list = list(c('Accuracy', 'Micro.F1'), c('Micro.Precision', 'Micro.Recall'), c('F1_5mC', 'F1_5C'))
for (measure.pair in measure.pair.list) {
  fig.34a.scatter.plot.performance(df, measure.pair, locations = locations.Singletons, bdir = outdir, scale = 0.65)
  fig.34a.scatter.plot.performance(df, measure.pair, locations = locations.Genome, bdir = outdir, scale = 0.75)
  #break
}

# Bar plot
source(here('src', 'rplot_func', 'utils_ggplot2_paper.R')) # Bar plot
perf.measure.list = c('Accuracy', 'ROC.AUC', 'Micro.F1', 'Macro.F1', 'Average.Precision', 'Recall_5mC')
for (perf.measure in perf.measure.list) {
  fig.34a.bar.plot.performance(df, perf.measure, bdir = outdir, locations = locations.Singletons, scale = 0.6)
  fig.34a.bar.plot.performance(df, perf.measure, bdir = outdir, locations = locations.Genome, scale = 0.75)
  #break
}

outfn = here('result', 'figure3a.work.RData')
save.image(file = outfn)
printf("save workspace env to %s", outfn)

