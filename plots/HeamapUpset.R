library(ComplexHeatmap)

lt = list(set1 = c("a", "b", "c"),
          set2 = c("b", "c", "d", "e"))
list_to_matrix(lt)

m = make_comb_mat(lt)
UpSet(m)

UpSet(m, set_order = c("a", "b", "c"), comb_order = order(comb_size(m)))

