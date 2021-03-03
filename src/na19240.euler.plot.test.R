rm(list = ls())

library(here)
library(eulerr)

source(here('src', 'plotutils4r', 'paper_utils.R'))

center.red = 3
reduce = 1
a = 56634168 / reduce
b = 53952133 / reduce
c = 56927603 / reduce
ab = 53567179 / reduce
ac = 56287734 / reduce
bc = 53790051 / reduce
abc = 53414153 / reduce

#comb = c("A" = 56634168, "B" = 53952133, "C" = 56927603,
#         "A&B" = 53567179, "A&C" = 56287734, "B&C" = 53790051,
#         "A&B&C" = 53414153)

comb = c("A" = a - ab - ac + abc, "B" = b - ab - bc + abc, "C" = c - ac - bc + abc,
         "A&B" = ab - abc, "A&C" = ac - abc, "B&C" = bc - abc,
         "A&B&C" = abc / center.red)

print(comb)

#fit1 <- euler(comb, input = 'union', shape = 'circle')
#fit1 <- euler(comb, shape = 'ellipse')
fit1 <- euler(comb, shape = 'circle')


print(fit1)

graphics.off()

p1 <- plot(fit1, quantities = list(type = c("counts"),
                                   font = 1, round = 2, cex = 0.6),
           fills = Top3.ToolColorPal,
           labels = identical(legend, FALSE),
           adjust_labels = TRUE,
           counts = TRUE,
           legend = list(labels = Top3.Tool.Order, cex = 0.5,
                         colors = Top3.ToolColorPal, alpha = 1, side = 'bottom', ncol = 3, nrow = 1, hgap = 0),
           alpha = 0.5,
           fill_opacity = 0.5,
           edges = list(lty = 1),
           pal = Top3.ToolColorPal
)

p1

out_dir = here('figures', 'venn-plot')
outfn = sprintf("%s/euller.plot.na19240.new.jpg", out_dir)
ggsave(p1, file = outfn, dpi = 600, width = 3, height = 2)