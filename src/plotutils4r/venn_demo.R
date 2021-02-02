library(VennDiagram)
library(ggplot2)
library(here)
library(eulerr)


source(here('src', 'rplot_func', 'paper_utils.R'))

outdir = here('figures')


fit1 <- euler(c("Day 1" = 69, "Day 2" = 109, "Day 3" = 152,
                "Day 1&Day 2" = 11, "Day 1&Day 3" = 18, "Day 2&Day 3" = 28,
                "Day 1&Day 2&Day 3" = 2), input = 'union', shape = 'ellipse')

col <- c("#66D2D6", "#FBC740", "#E56997")

plot(fit1,
     #quantities = TRUE,
     fills = list(fill = col, alpha = 0.6),
     edges = "black",
     quantities = list(type = c("counts", 'percent'), cex = 1.5),
     #counts=list(cex=5),
     main = list(label = "test", cex = 2),
     legend = list(side = "bottom", nrow = 1, ncol = 3, vgap = 5, hgap=0.5),
     adjust_labels = TRUE)


# Venn diagram plot
for (fn in venn_flist) {
  infn = here('result', fn)
  dt <- read.table(infn)

  base_infn = basename(infn)
  outfn = sprintf("%s/%s.jpg", outdir, base_infn)

  fig.6a.venn.plot.set5(dt$V1, outfn)
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


library(eulerr)

# Input in the form of a named numeric vector
fit1 <- euler(c("A" = 25, "B" = 5, "C" = 5,
                "A&B" = 5, "A&C" = 5, "B&C" = 3,
                "A&B&C" = 3))


ret = dt$V1
fit1 <- euler(c("A" = ret[1], "B" = ret[2], "C" = ret[3],
                "A&B" = ret[4], "A&C" = ret[5], "B&C" = ret[6],
                "A&B&C" = ret[7]), input = 'union', shape = 'circle')


gp <- plot(fit1,
           quantities = list(type = c("counts", "percent"),
                             font = 1, round = 2, cex = 1),
           labels = NULL,
           pal = Top3.ToolColorPal,
           fills = Top3.ToolColorPal,
           alpha = 0.5,
           fill_opacity = 0.5,
           edges = list(lty = 1),
           main = 'Title',
           counts = TRUE,
           legend = list(labels = Top3.Tool.Order, cex = 1.5)
)

ggsave(gp, file = 'test-euller.jpg', dpi = 600)

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
