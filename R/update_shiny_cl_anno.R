options(stringsAsFactors = FALSE)


# Load data
# feather.path <- "//allen/programs/celltypes/workgroups/hct/CT_clustering/mouse_cell_vs_nuc/shinydb/20170818_VISp_L5_cell_1/"
# anno <- as.data.frame(read_feather(paste0(feather.path, "anno.feather")))
# 
# cl.info <- unique(anno[, c("cluster_label", "cluster_type_label")])
# cl.info$cluster_color <- NA
# 
# for (type1 in unique(cl.info$cluster_type_label)) {
#   if (type1 %in% c("donor", "outlier")) {
#     cl.colors <- grey.colors(length(cl.info$cluster_color[cl.info$cluster_type == type1]))
#   } else if (type1 == "glia") {
#     cl.colors <- MakeColors(1:length(cl.info$cluster_color[cl.info$cluster_type == type1]),
#                             h.range = c(0, 360), c.range = c(0, 0.2), l.range = c(0, 1.5))
#   } else if (type1 == "inh") {
#     cl.colors1 <- MakeColors(1:round(length(cl.info$cluster_color[cl.info$cluster_type == type1]) / 2, 0),
#                              h.range = c(0, 80), c.range = c(0, 1.5), l.range = c(0.5, 1.5))
#     cl.colors2 <- MakeColors(1:(length(cl.info$cluster_color[cl.info$cluster_type == type1]) - length(cl.colors1)), c(280, 360), c(0, 1.5), c(0.5, 1.5))
#     cl.colors <- sample(c(cl.colors1, cl.colors2), length(c(cl.colors1, cl.colors2)), replace = FALSE)
#   } else if (type1 == "exc") {
#     cl.colors <- MakeColors(1:length(cl.info$cluster_color[cl.info$cluster_type == type1]),
#                             c(100, 260), c(0, 1.5), c(0.5, 1.5))
#   }
#   cl.info$cluster_color[cl.info$cluster_type == type1] <- cl.colors
# }
# 
# 
# for (cl1 in unique(anno$cluster_label)) {
#   anno$cluster_color[anno$cluster_label == cl1] <- cl.info$cluster_color[cl.info$cluster_label == cl1]
# }
# 
# write_feather(anno, paste0(feather.path, "anno.feather"))
