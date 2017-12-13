library(pheatmap)
library(RColorBrewer)
source("//allen/programs/celltypes/workgroups/rnaseqanalysis/Script_Repository/Trygve/src/fOrderConfusion.R")

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_nuclei/cl.rda",
     verbose = TRUE)
cl.z <- data.frame(cl, cellid = names(cl))

cl.df <- read.csv("//allen/programs/celltypes/workgroups/hct/CT_clustering/mouse_cell_vs_nuc/cache/20171211_VISp_nuc/1/cl_merge.csv", row.names = 1)
cl.df$cellid <- sapply(rownames(cl.df), function(x) strsplit(x, "~")[[1]][1])

cl.df.all <- merge(cl.z, cl.df, by = "cellid")

cc1 <- OrderConfusion(cl.df.all$cl, cl.df.all$x)

hm.colors <- colorRampPalette(c("white", brewer.pal(9, "YlOrRd")))(100)
pheatmap(cc1, cluster_rows = FALSE, cluster_cols = FALSE, color = hm.colors)
