library(ComplexHeatmap)

lt = list(set1 = c("a", "b", "c"),
          set2 = c("b", "c", "d", "e"))
list_to_matrix(lt)

m = make_comb_mat(lt)

up1 = UpSet(m)
png(file="filename_heatmap1.png")

# UpSet(m, set_order = c("a", "b", "c"), comb_order = order(comb_size(m)))

draw(up1)

dev.off()



movies = read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
                  header = TRUE, sep = ";")
head(movies) # `make_c

m = make_comb_mat(movies, top_n_sets = 6)
m



