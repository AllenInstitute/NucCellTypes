library(feather)
library(pheatmap)



all.paths <- c("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20170905_exons_introns/", 
               "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20170405/",
               "//allen/programs/celltypes/workgroups/hct/CT_clustering/mouse_cell_vs_nuc/shinydb/20170818_VISp_L5_nuc_1/",
               "//allen/programs/celltypes/workgroups/hct/CT_clustering/mouse_cell_vs_nuc/shinydb/20170818_VISp_L5_cell_1/",
               "//allen/programs/celltypes/workgroups/hct/CT_clustering/mouse_cell_vs_nuc/shinydb/20170818_VISp_L5_snap25_1/",
               "//allen/programs/celltypes/workgroups/hct/CT_clustering/mouse_cell_vs_nuc/shinydb/20170818_VISp_L5_nuc_exon_1/")
names(all.paths) <- c("Mouse_V1_ALM", "Mouse_V1_ALM_0405", 
                      "Mouse_L5_Nuc", "Mouse_L5_Cell", 
                      "Mouse_L5_Snap25", "Mouse_L5_Nuc_Ex")
paths <- all.paths

dend <- list()
anno <- list()
anno.cl <- list()
prop.all <- list()
for (dat1 in names(paths)) {
  if (file.exists(paste0(paths[[dat1]], "dend.RData"))) {
    dend[[dat1]] <- readRDS(paste0(paths[[dat1]], "dend.RData"))
  }
  anno[[dat1]] <- as.data.frame(read_feather(paste0(paths[[dat1]], "anno.feather")))
  colnames(anno[[dat1]])[colnames(anno[[dat1]]) == "final_label"] <- "cluster_label"  # KLUDGE - Tasic
  anno[[dat1]]$cluster_label <- sub("_outlier", "", anno[[dat1]]$cluster_label)  # KLUDGE - fix MTG labels
  prop1 <- as.data.frame(t(read_feather(paste0(paths[[dat1]], "prop.feather"))))
  rownames(prop1) <- toupper(rownames(prop1))
  
  if (dat1 == "Mouse_V1_ALM") {
    cl.lab <- unique(anno[[dat1]][, c("cluster_id", "cluster_label")])
    cl.lab <- cl.lab[order(cl.lab$cluster_id), ]
    colnames(prop1) <- as.character(cl.lab$cluster_label)
  } else {
    colnames(prop1) <- unique(anno[[dat1]]$cluster_label)
  }
  keep.cl <- unlist(dendrapply(dend[[dat1]], function(x) if (is.leaf(x)) attr(x, "label")))
  # keep.cl <- unique(anno[[dat1]][anno[[dat1]]$cluster_type_label %in% c("inh", "exc", "glia"), c("cluster_label")])
  anno.cl1 <- data.frame(cluster_label = unique(anno[[dat1]][, c("cluster_label")]))
  
  cl.size <- table(anno[[dat1]]$cluster_label)
  anno.cl1$size <- cl.size[match(names(cl.size), anno.cl1$cluster_label)]
  anno.cl[[dat1]] <- anno.cl1[match(keep.cl, anno.cl1$cluster_label), ]
  anno[[dat1]] <- droplevels(subset(anno[[dat1]], cluster_label %in% keep.cl))
  prop1.subset <- prop1[, match(keep.cl, colnames(prop1))]
  prop.all[[dat1]] <- prop1.subset
}


# Cre line proportions of cells matched to nuclei
cre.cnt <- table(anno[["Mouse_L5_Cell"]]$cluster_label, 
                 anno[["Mouse_L5_Cell"]]$cre_line_label)
write.csv(cre.cnt, file = "output/cre_cluster.csv")
cre.prop <- sweep(cre.cnt, 1, rowSums(cre.cnt), "/")
pheatmap(cre.prop, cluster_rows = FALSE)


# Layer proportions of clusters
layer.cnt <- table(anno[["Mouse_V1_ALM"]]$cluster_label, 
                 anno[["Mouse_V1_ALM"]]$layer_label)
layer.prop <- sweep(layer.cnt, 1, rowSums(layer.cnt), "/")
pheatmap(layer.prop, cluster_cols = FALSE)

layer.prop[grep("Foxp2_1", rownames(layer.prop)), ]

# Missing nuclear cluster
layer.cnt <- with(droplevels(subset(anno[["Mouse_V1_ALM"]], cluster_label == "L6a Foxp2_1")),
                  table(cre_label, layer_label))
layer.prop <- sweep(layer.cnt, 1, rowSums(layer.cnt), "/")
print(layer.prop)


# Layer proportions of Cre line sampling
layer.cnt <- table(anno[["Mouse_V1_ALM"]]$cre_label, 
                   anno[["Mouse_V1_ALM"]]$layer_label)
layer.prop <- sweep(layer.cnt, 1, rowSums(layer.cnt), "/")
pheatmap(layer.prop, cluster_rows = FALSE, cluster_cols = FALSE)

cre.prop <- sweep(layer.cnt, 2, colSums(layer.cnt), "/")
pheatmap(cre.prop, cluster_rows = FALSE, cluster_cols = FALSE)
