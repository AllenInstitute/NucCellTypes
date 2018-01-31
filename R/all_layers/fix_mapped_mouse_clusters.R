library(feather)
source("//allen/programs/celltypes/workgroups/rnaseqanalysis/Script_Repository/Trygve/src/celltypes/fMaptoV1Tree.R")


root.path <- "//allen/programs/celltypes/workgroups"
shiny.path <- paste0(root.path, "/hct/CT_clustering/mouse_cell_vs_nuc/shinydb/20171211_VISp_nuc_1/")
anno <- as.data.frame(read_feather(paste0(shiny.path, "anno.feather")))
data <- as.data.frame(read_feather(paste0(shiny.path, "data_t.feather")))
row.names(data) <- data$gene
data <- data[, -grep("gene", colnames(data))]

#### Map to mouse ####
max.leaf <- apply(data, 2, 
                  function(x) MaptoV1Tree(expr = x, 
                                          genes = rownames(data),
                                          root.path = root.path)$leaf.max)
max.leaf.class <- apply(data, 2, 
                        function(x) MaptoV1Tree(expr = x, 
                                                genes = rownames(data),
                                                root.path = root.path)$leaf.max.class)




#### Define mapped V1 cluster colors/order ####
type.dat <- read.csv(file = paste0(root.path, "/humancelltypes/Analyses/CellTypes/V1_mouse/lib/v1_cluster_metadata_Tasic.csv"))

order.map <- match(max.leaf, make.names(type.dat$Tree_leaf))
order.map[is.na(order.map)] <- match(max.leaf, make.names(type.dat$Tasic_et_al_2016_label))[is.na(order.map)]
order.map[is.na(order.map)] <- match(max.leaf, make.names(type.dat$vignette_label))[is.na(order.map)]

map.col <- type.dat$vignette_color[order.map]
map.cl.order <- type.dat$cluster_order[order.map]

anno$cellmap_id <- map.cl.order
anno$cellmap_label <- max.leaf
anno$cellmap_color <- map.col
anno$max.leaf_id <- map.cl.order
anno$max.leaf_label <- max.leaf
anno$max.leaf_color <- map.col
anno$max.leaf.class_id <- as.numeric(as.factor(max.leaf.class))
anno$max.leaf.class_label <- max.leaf.class
anno$max.leaf.class_color <- c("blue", "orange", "grey")[anno$max.leaf.class_id]

write_feather(anno, paste0(shiny.path, "anno.feather"))
