rm(list = ls())

library(here)
source(here('src', 'plotutils4r', 'paper_utils.R'))


# Load data and sort string orders
source(here('src', 'plotutils4r', 'paper_utils.R'))
infn = here('result', 'performance-results.csv')
df <- load.performance.data(infn)

## Figure 3, 4 a:Line plot
source(here('src', 'plotutils4r', 'paper_utils.R'))
out_dir = here('figures', 'line-plot')
dir.create(out_dir, showWarnings = FALSE)

for (perf.measure in measure.list) {
  fig.34a.line.plot.performance(df, perf.measure, bdir = out_dir, locations = locations.Singletons)
  fig.34a.line.plot.performance(df, perf.measure, bdir = out_dir, locations = locations.Genome)
  #break
}


## Figure 3, 4 b:Box plot
source(here('src', 'plotutils4r', 'paper_utils.R'))
out_dir = here('figures', 'box-plot')
dir.create(out_dir, showWarnings = FALSE)
for (perf.measure in measure.list) {
  fig.34b.box.location.performance(df, perf.measure, bdir = out_dir, locations = locations.Singletons)
  fig.34b.box.location.performance(df, perf.measure, bdir = out_dir, locations = locations.Genome)
  #break
}

## Figure 5 c: Resource summary
source(here('src', 'plotutils4r', 'paper_utils.R'))
outdir = here('figures')
fig.5c.running.resource.bar.plot(outdir)


## Figure 6 ab: Set venn and euller plot
source(here('src', 'plotutils4r', 'paper_utils.R'))
data_dir = here('result', 'venn-data')
out_dir = here('figures', 'venn-plot')
dir.create(out_dir, showWarnings = FALSE)
pattern.str = 'venn.data.*.dat'
pattern.str = 'venn.data.NA19240.*.dat'
for (venfn in list.files(data_dir, pattern = pattern.str)) {
  infn = here('result', 'venn-data', venfn)
  dt <- read.table(infn)
  base_infn = basename(infn)

  if (length(dt$V1) == 31) {
    #outfn = sprintf("%s/venn.plot.%s.jpg", out_dir, base_infn)
    #fig.6a.venn.plot.set5(dt$V1, outfn)
  }else if (length(dt$V1) == 7) {
    outfn2 = sprintf("%s/euller.plot.%s.jpg", out_dir, base_infn)
    fig.6b.euller.plot.set3(dt$V1, outfn2)
    #break
  }
  #break
}


## Figure 5d: Violin plot of corr COE
source(here('src', 'plotutils4r', 'paper_utils.R'))
outdir = here('figures')
fig.5d.violin.corr.performance(df, outdir)


## Figure 3, 4 a:Bar plot
source(here('src', 'plotutils4r', 'paper_utils.R'))
for (measure.pair in measure.pair.list) {
  for (dsname in Dataset.Order) {
    fig.34a.bar.dataset.location.performance(df, dsname, measure.pair, outdir)
  }
  #break
}


quit()

outfn = here('result', 'figure3a.work.RData')
save.image(file = outfn)
printf("save workspace env to %s", outfn)

