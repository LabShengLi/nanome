rm(list = ls())
library(here)

# Figure 5b: Correlation bar for each region
source(here('src', 'plotutils4r', 'paper_utils.R'))
library(readxl)
library(tidyverse)
out_dir = here('figures')
infn = here('result', 'All.corrdata.coe.pvalue.each.regions.xlsx')
df = read_excel(infn)

df$dsname <- factor(df$dsname, levels = Dataset.Order)
df$Tool <- factor(df$Tool, levels = Tool.Order[1:5])
df$Location <- factor(df$Location, levels = c(locations.Singletons, locations.Genome))
df <- df %>%
  drop_na(Tool, Location) %>%
  filter(dsname == 'NA19240')

p1 <- ggplot(data = df, aes(x = Tool, y = COE, fill = Tool)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~Location, ncol = 5) +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(strip.text.x = element_text(size = 10)) +
  scale_fill_manual(values = ToolColorPal) +
  ylab("Pearson correlation coefficient") +
  theme(legend.title = element_blank()) +
  ylim(min(df$COE), 1)

p1

outfn = sprintf("%s/bar.plot.coe.in.each.regions.jpg", out_dir)
ggsave(p1, filename = outfn, width = 6.5, height = 5, dpi = 600)
printf("save to %s\n", outfn)


# Load data and sort string orders
source(here('src', 'plotutils4r', 'paper_utils.R'))
infn = here('result', 'performance-results-cut5.csv')
df <- load.performance.data(infn)

## Figure S for 3, 4: Line plot
source(here('src', 'plotutils4r', 'paper_utils.R'))
out_dir = here('figures', 'line-plot')
dir.create(out_dir, showWarnings = FALSE)

for (perf.measure in measure.list) {
  fig.s34a.line.plot.performance(df, perf.measure, bdir = out_dir, locations = locations.Singletons)
  fig.s34a.line.plot.performance(df, perf.measure, bdir = out_dir, locations = locations.Genome)
  #break
}


## Figure 3, 4 a:Box plot
source(here('src', 'plotutils4r', 'paper_utils.R'))
out_dir = here('figures', 'box-plot')
dir.create(out_dir, showWarnings = FALSE)
for (perf.measure in measure.list) {
  fig.34a.box.location.performance(df, perf.measure, bdir = out_dir, locations = locations.Singletons)
  fig.34a.box.location.performance(df, perf.measure, bdir = out_dir, locations = locations.Genome)
  #break
}

## Figure 5 c: Resource summary
source(here('src', 'plotutils4r', 'paper_utils.R'))
outdir = here('figures')
fig.5c.running.resource.bar.plot(outdir)


## Figure 5d: Violin plot of corr COE
source(here('src', 'plotutils4r', 'paper_utils.R'))
outdir = here('figures')
fig.5d.violin.corr.performance(df, outdir)


## Figure 6 ab: Set venn and euller plot
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
    fig.6a.venn.plot.set5(dt$V1, outfn)
  }else if (length(dt$V1) == 7) {
    outfn2 = sprintf("%s/euller.plot.%s.jpg", out_dir, base_infn)
    fig.6b.euller.plot.set3(dt$V1, outfn2)
    #break
  }
  #break
}


## Figure 6 c: bar plot of number of sites for each tool
source(here('src', 'plotutils4r', 'paper_utils.R'))
fig.6c.bar.plot.tools()


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

