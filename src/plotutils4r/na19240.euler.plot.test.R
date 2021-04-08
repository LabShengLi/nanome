rm(list = ls())

library(here)
library(eulerr)

source(here('src', 'plotutils4r', 'paper_utils.R'))
# a-DeepSignal b-Nanopolish c-Megalodon
center.red = 1
reduce = 1

# modified data for clear regions
a = 57634168 / reduce
b = 53576718 / reduce
c = 57927603 / reduce
ab = 52825344 / reduce
ac = 56287734 / reduce
bc = 53033412 / reduce
abc = 52687665 / reduce

# original data
#a = 56634168 / reduce
#b = 53176718 / reduce
#c = 56927603 / reduce
#ab = 52825344 / reduce
#ac = 56287734 / reduce
#bc = 53033412 / reduce
#abc = 52687665 / reduce

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

# Original data results
#print(comb)
#       A        B        C      A&B      A&C      B&C    A&B&C
#  208755     5627   294122   137679  3600069   345747 52687665

#   original      fitted  residuals regionError
#A       208755   212373.31  -3618.309       0.000
#B         5627    24200.31 -18573.310       0.000
#C       294122   303526.66  -9404.663       0.000
#A&B     137679        0.00 137679.000       0.002
#A&C    3600069  3599537.15    531.854       0.000
#B&C     345747   338775.07   6971.934       0.000
#A&B&C 52687665 52687625.71     39.291       0.002


graphics.off()

#quantities = list(type = c("counts"),
#                                   font = 1, round = 2, cex = 0.6),
p1 <- plot(fit1,
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
ggsave(p1, file = outfn, dpi = 600, width = 2, height = 2)