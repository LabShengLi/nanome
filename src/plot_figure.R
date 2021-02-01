rm(list = ls())

library(here)
here()

getwd()
#wdir = here('src')  #"/Users/liuya/PycharmProjects/nano-compare/src/rplot_func"
#setwd(wdir)

source(here('src', 'plotutils4r', 'paper_utils.R'))
outdir = here('figures')

fig5c.running.resource.bar.plot(outdir)

# Load data and sort string orders
infn = here('result', 'performance-results.csv')
df <- load.performance.data(infn)

## Test
source(here('src', 'plotutils4r', 'paper_utils.R'))

for (corr_col in Corr.Perf.List) {
  fig.34a.violin.corr.performance(df, corr_col, outdir, scale = 1)
}


## Set venn and euller plot
source(here('src', 'plotutils4r', 'paper_utils.R'))
data_dir = here('result', 'venn-data')
out_dir = here('figures', 'venn-plot')
dir.create(out_dir, showWarnings = FALSE)
for (venfn in list.files(data_dir, 'venn.data.*.dat')) {
  infn = here('result', 'venn-data', venfn)
  dt <- read.table(infn)
  base_infn = basename(infn)

  if (length(dt$V1) == 31) {
    outfn = sprintf("%s/venn.plot.%s.jpg", out_dir, base_infn)
    fig.34c.venn.plot.set5(dt$V1, outfn)
  }else if (length(dt$V1) == 7) {
    #outfn1 = sprintf("%s/venn.plot.%s.jpg", out_dir, base_infn)
    #fig.34c.venn.plot.set3(dt$V1, outfn1)
    outfn2 = sprintf("%s/euller.plot.%s.jpg", out_dir, base_infn)
    fig.34c.euller.plot.set3(dt$V1, outfn2)
    #break
  }
  #break
}


## Box plot
for (perf.measure in measure.list) {
  fig.34a.box.location.performance(df, perf.measure, bdir = outdir, locations = locations.Singletons)
  fig.34a.box.location.performance(df, perf.measure, bdir = outdir, locations = locations.Genome)
  #break
}


## Figure 3, 4 a:Bar plot
source(here('src', 'plotutils4r', 'paper_utils.R'))

for (measure.pair in measure.pair.list) {
  for (dsname in Dataset.Order) {
    fig.34a.bar.dataset.location.performance(df, dsname, measure.pair, outdir)
  }
  #break
}

## Figure 3, 4 a:Line plot
source(here('src', 'plotutils4r', 'paper_utils.R'))
out_dir = here('figures', 'line-plot-CPG')
dir.create(out_dir, showWarnings = FALSE)

for (perf.measure in measure.list) {
  fig.34a.line.plot.performance(df, perf.measure, bdir = out_dir, locations = locations.Singletons)
  fig.34a.line.plot.performance(df, perf.measure, bdir = out_dir, locations = locations.Genome)
  #break
}


quit()

outfn = here('result', 'figure3a.work.RData')
save.image(file = outfn)
printf("save workspace env to %s", outfn)

