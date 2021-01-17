library(VennDiagram)
library(ggplot2)

# NOT RUN {
# Reference five-set diagram
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
  category = c("DeepSignal", "Tombo", "Nanopolish", "DeepMod", "Megalodon"),
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.cex = 2,
  margin = 0.05,
  cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
          1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
  ind = TRUE
);

# Writing to file
#tiff(filename = "Quintuple_Venn_diagram.tiff", compression = "lzw");
ggsave(venn.plot, filename = "Quintuple_Venn_diagram.tiff")
#grid.draw(venn.plot);
#dev.off();
# }