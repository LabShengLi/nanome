rm(list = ls())
library(here)

library(UpSetR)


source(here('src', 'plotutils4r', 'paper_utils.R'))
data_dir = here('result', 'venn-data')
out_dir = here('figures', 'venn-plot')
dir.create(out_dir, showWarnings = FALSE)

pattern.str = 'venn.data.*.dat'
# pattern.str = "venn.data.NA19240_RRBS_2Reps.NA19240.five.tools.cov3.dat"
# pattern.str = "venn.data.APL_RRBS.APL.five.tools.cov3.dat"
pattern.str = "venn.data.HL60_RRBS_2Reps.HL60.five.tools.cov3.dat"
print(list.files(data_dir, pattern = pattern.str))

for (venfn in list.files(data_dir, pattern = pattern.str)) {
    infn = here('result', 'venn-data', venfn)
    print(infn)
    dt <- read.table(infn)
    base_infn = basename(infn)
    
    if (length(dt$V1) == 7) next
    
    # if (length(dt$V1) == 31) {
    outfn = sprintf("%s/upset.plot.%s.pdf", out_dir, base_infn)
    fit1 <- fig.5newA.upset.plot.set5(dt$V1, outfn)
    
    # graphics.off()
    # dev.off()
    print("start upset plotting")
    
    p1 <- upset(fromExpression(fit1$original.values), scale.intersections = "identity",
                sets.x.label = "Total CpGs", mainbar.y.label = "Intersection CpGs",
                show.numbers = 'yes', text.scale = 2, number.angles = 30, nintersects = 31, keep.order = TRUE,
                sets.bar.color = "grey", set_size.show = TRUE, set_size.angles = 0, set_size.numbers_size = 6,
                order.by = "freq")
    
    pdf(file = outfn, onefile = FALSE) # or other device
    p1
    dev.off()
    printf('save to %s', outfn)
    
    
    pdf(file = outfn, onefile = FALSE) # or other device
    p1
    print("Start saving")
    dev.off()
    printf('save to %s', outfn)
    
    break
    # }else
}


# Figure S1: pie plot of singleton and nonsingleton sites
source(here('src', 'plotutils4r', 'paper_utils.R'))
fig.s1.pie.plot.singletons.nonsingletons.raw.fast5()

# Export Tabls S3 by order
source(here('src', 'plotutils4r', 'paper_utils.R'))
export.table.s3.xlsx()

# Export Tabls S4 by order
source(here('src', 'plotutils4r', 'paper_utils.R'))
export.table.s4.xlsx()

## Figure S for 3, 4: Line plot
source(here('src', 'plotutils4r', 'paper_utils.R'))
out_dir = here('figures', 'line-plot')
dir.create(out_dir, showWarnings = FALSE)

for (perf.measure in measure.list) {
    fig.s34a.line.plot.performance(perf.measure, bdir = out_dir, locations = locations.Singletons)
    fig.s34a.line.plot.performance(perf.measure, bdir = out_dir, locations = locations.Regions)
    #break
}

## Figure 3, 4 a:Box plot
source(here('src', 'plotutils4r', 'paper_utils.R'))
out_dir = here('figures', 'box-plot')
dir.create(out_dir, showWarnings = FALSE)
for (perf.measure in measure.list) {
    fig.34a.box.location.performance(perf.measure, bdir = out_dir, locations = locations.Singletons)
    fig.34a.box.location.performance(perf.measure, bdir = out_dir, locations = locations.Regions)
    #break
}

# Figure 5b: Correlation COE bar for each region
source(here('src', 'plotutils4r', 'paper_utils.R'))
fig.5b.sorted.bar.plot.coe.in.each.region()


## Figure 7b1: Resource summary
source(here('src', 'plotutils4r', 'paper_utils.R'))
fig.7b1.running.resource.bar.plot()


## Figure 7b2: Benchmarking
source(here('src', 'plotutils4r', 'paper_utils.R'))
fig.7b2.running.resource.bar.plot()


## Figure 5 cd: Set venn and euller plot
source(here('src', 'plotutils4r', 'paper_utils.R'))
data_dir = here('result', 'venn-data')
out_dir = here('figures', 'venn-plot')
dir.create(out_dir, showWarnings = FALSE)
pattern.str = 'venn.data.*.dat'
#pattern.str = 'venn.data.NA19240.*.dat'
for (venfn in list.files(data_dir, pattern = pattern.str)) {
    infn = here('result', 'venn-data', venfn)
    dt <- read.table(infn)
    base_infn = basename(infn)
    
    if (length(dt$V1) == 31) {
        outfn = sprintf("%s/venn.plot.%s.jpg", out_dir, base_infn)
        fig.5cd.venn.plot.set5(dt$V1, outfn)
    }else if (length(dt$V1) == 7) {
        outfn2 = sprintf("%s/euller.plot.%s.jpg", out_dir, base_infn)
        fig.5cd.euller.plot.set3(dt$V1, outfn2)
        #break
    }
    #break
}


## Figure 6A: bar plot of number of sites for each tool
source(here('src', 'plotutils4r', 'paper_utils.R'))
fig.6a.bar.plot.tools.sites.all.datasets()


quit()

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


## Figure 5d: Violin plot of corr COE
source(here('src', 'plotutils4r', 'paper_utils.R'))
outdir = here('figures')
fig.5d.violin.corr.performance(df, outdir)

