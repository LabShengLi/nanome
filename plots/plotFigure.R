rm(list = ls())
library(here)
plotsDir = "plots"
resultDir = "result"
figuresDir = "figures"
utilFileName = "plotUtils.R"

########################################
########################################
# Figure S1: pie plot of singleton and nonsingleton sites
source(here(plotsDir, utilFileName))
fig.s1.pie.plot.singletons.nonsingletons.raw.fast5()

########################################
########################################
# Export Tabls S3 by order
source(here(plotsDir, utilFileName))
export.table.s3.xlsx()

########################################
########################################
# Export Tabls S5 by order
source(here(plotsDir, utilFileName))
export.table.s5.xlsx()

########################################
########################################
## Figure S3, S4: Line plot
source(here(plotsDir, utilFileName))
out_dir = here(figuresDir, 'line-plot')
dir.create(out_dir, showWarnings = FALSE)

for (perf.measure in measure.list) {
    fig.s34a.line.plot.performance(perf.measure, out_dir = out_dir, locations = locations.Singletons)
    fig.s34a.line.plot.performance(perf.measure, out_dir = out_dir, locations = locations.Regions)
    #break
}

########################################
########################################
## Figure 3a,4a,S34b:Box plot
source(here(plotsDir, utilFileName))
out_dir = here(figuresDir, 'box-plot')
dir.create(out_dir, showWarnings = FALSE)

for (perf.measure in measure.list) {
    fig.34a.box.location.performance(perf.measure, out_dir = out_dir, locations = locations.Singletons)
    fig.34a.box.location.performance(perf.measure, out_dir = out_dir, locations = locations.Regions)
    #break
}


########################################
########################################
## Figure 6: Venn and euller plot
source(here(plotsDir, utilFileName))
data_dir = here(resultDir, 'venn-data')
out_dir = here(figuresDir, 'venn-plot')
dir.create(out_dir, showWarnings = FALSE)

pattern.str = 'venn.data.*.dat'
flist = list.files(data_dir, pattern = glob2rx(pattern.str))
print(flist)
for (venfn in flist) {
    infn = here(resultDir, 'venn-data', venfn)
    dt <- read.table(infn)
    base_infn = basename(infn)
    if (length(dt$V1) == 31) { #next # No Venn plot needed
        outfn = sprintf("%s/fig.6ab.venn.plot.%s.jpg", out_dir, base_infn)
        fig.6ab.venn.plot.set5(dt$V1, outfn)
    }else if (length(dt$V1) == 7) { ##next
        outfn2 = sprintf("%s/fig.6ab.euller.plot.%s.jpg", out_dir, base_infn)
        fig.6ab.euller.plot.set3(dt$V1, outfn2)
        #break
    }
    #break
}


########################################
########################################
## Figure S6: sorted bar of COE accross regions
source(here(plotsDir, utilFileName))
fig.s6.sorted.bar.plot.coe.in.each.region()

########################################
########################################
## Figure 7, Table S7: Resource summary
source(here(plotsDir, utilFileName))
fig.7b1.running.resource.bar.plot()

print("plotFigure DONE")

quit()

########################################
########################################
## Figure 7b2: Benchmarking
source(here('../src', 'plotutils4r', 'paper_utils.R'))
fig.7b2.benchmarking.running.resource.bar.plot()

## Figure S6b: bar plot of number of sites, deprecated
source(here(plotsDir, utilFileName))
fig.6a.bar.plot.tools.sites.all.datasets()

## Figure 3, 4 a:Bar plot
source(here('../src', 'plotutils4r', 'paper_utils.R'))
for (measure.pair in measure.pair.list) {
    for (dsname in Dataset.Order) {
        fig.34a.bar.dataset.location.performance(df, dsname, measure.pair, outdir)
    }
    #break
}


quit()

outfn = here('../result', 'figure3a.work.RData')
save.image(file = outfn)
printf("save workspace env to %s", outfn)


## Figure 5d: Violin plot of corr COE
source(here('../src', 'plotutils4r', 'paper_utils.R'))
outdir = here('../figures')
fig.5d.violin.corr.performance(df, outdir)
