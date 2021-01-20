library(VennDiagram)
library(ggplot2)
library(here)


source(here('src', 'rplot_func', 'utils_ggplot2_paper.R'))

outdir = here('figures')

# Venn diagram plot
for (fn in venn_flist) {
  infn = here('result', fn)
  dt <- read.table(infn)

  base_infn = basename(infn)
  outfn = sprintf("%s/%s.jpg", outdir, base_infn)

  fig.34c.venn.plot(dt$V1, outfn)
}


graphics.off()


venn.plot <- draw.quintuple.venn(
  area1 = 301,
  area2 = 321,
  area3 = 311,
  area4 = 321,
  area5 = 301,
  n12 = 188,
  n13 = 191,
  n14 = 184,
  n15 = 177,
  n23 = 194,
  n24 = 197,
  n25 = 190,
  n34 = 190,
  n35 = 173,
  n45 = 186,
  n123 = 112,
  n124 = 108,
  n125 = 108,
  n134 = 111,
  n135 = 104,
  n145 = 104,
  n234 = 111,
  n235 = 107,
  n245 = 110,
  n345 = 100,
  n1234 = 61,
  n1235 = 60,
  n1245 = 59,
  n1345 = 58,
  n2345 = 57,
  n12345 = 31,
  category = Tool.Order,
  fill = ToolColorPal[1:5],
  cat.col = ToolColorPal[1:5],
  cat.pos = c(0, 2, 8, 5, 11),
  cat.dist = c(0.22, 0.22, -0.2, -0.2, 0.22),
  cat.cex = 2,
  margin = 0.05,
  cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
          1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
  ind = TRUE
);


outfn = sprintf("%s/ven-demo.pdf", bdir)

#tiff(filename = outfn, compression = "lzw");
#grid.draw(venn.plot);
#dev.off();

ggsave(venn.plot, file = outfn, dpi = 600)

printf('save to %s', outfn)



venn.plot <- draw.triple.venn(
	area1 = 65,
	area2 = 75,
	area3 = 85,
	n12 = 35,
	n23 = 15,
	n13 = 25,
	n123 = 5,
	category = c("First", "Second", "Third"),
	fill = c("blue", "red", "green"),
	lty = "blank",
	cex = 2,
	cat.cex = 2,
	cat.col = c("blue", "red", "green")
	);
