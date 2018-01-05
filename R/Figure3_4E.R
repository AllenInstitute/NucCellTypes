library(feather)
library(plyr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)

# Load nuc/cell matching
nuc.ids <- scan("data/nuc_ids.txt", "character")
cell.ids <- scan("data/cell_ids.txt", "character")

# Load cc heatmaps
dirs <- c("nuc", "nuc_varE_clIE", "nuc_varIE_clE", "nuc_exon", 
          "cell", "cell_varE_clIE", "cell_varIE_clE", "cell_exon")
sets <- c("nuc_IE_IE", "nuc_E_IE", "nuc_IE_E", "nuc_E_E", 
          "cell_IE_IE", "cell_E_IE", "cell_IE_E", "cell_E_E")
names(dirs) <- sets
paths <- paste0("data/20170818_VISp_L5_", dirs)
names(paths) <- sets


cc.list <- list()
anno.list <- list()
for (set1 in sets) {
  anno.list[[set1]] <- as.data.frame(read_feather(paste0(paths[set1], "/anno.feather")))
  cc1 <- read.csv(file = paste0(paths[set1], "/cl.cons.csv.gz"),
                  row.names = 1)
  cc.ids <- sapply(row.names(cc1), function(x) strsplit(x, "~")[[1]][1])
  cc.order <- match(anno.list[[set1]]$sample_id, cc.ids)
  cc.list[[set1]] <- cc1[cc.order, cc.order]
}



#### Figure 3A,B - Within/between co-clustering - all comparisons ####
hm.colors <- colorRampPalette(c("white", brewer.pal(9, "YlOrRd")))(100)
for (sample.type1 in c("nuc", "cell")) {
  for (set1 in sets[grep(sample.type1, sets)]) {
    if (set1 == sets[grep(sample.type1, sets)][1]) {
      ph.ie <- pheatmap(cc.list[[set1]], silent = TRUE)
      ie.order <- ph.ie$tree_row$order
      ph.exon <- pheatmap(cc.list[[sets[grep(sample.type1, sets)][4]]], silent = TRUE)
      exon.order <- ph.exon$tree_row$order
      dend1 <- as.dendrogram(ph.ie$tree_row)
      # Match exon only heatmap
      dend1.reordered <- reorder(dend1, match(ph.exon$tree_row$labels[exon.order], 
                                              ph.ie$tree_row$labels[ie.order]),
                                 FUN = "mean")
      hm.order2 <- as.hclust(dend1.reordered)$order
    }
    png(file = paste0("output/", set1, ".png"), 
        width = 1024, height = 1024, res = 300)
    pheatmap(cc.list[[set1]][hm.order2, hm.order2], 
             cluster_rows = FALSE, cluster_cols = FALSE, main = set1, 
             color = hm.colors, show_rownames = FALSE, show_colnames = FALSE)
    dev.off()    
  }
}



# Within/between co-clustering stats
cc.stats <- data.frame()
for (sample.type1 in c("nuc", "cell")) {
  for (set1 in sets[grep(sample.type1, sets)]) {
    cl.lab <- anno.list[[sets[grep(sample.type1, sets)][1]]]$cluster_label
    cc.mean1 <- apply(cc.list[[set1]], 1, function(x) tapply(x, cl.lab, mean))
    cc.mean <- apply(cc.mean1, 1, function(x) tapply(x, cl.lab, mean))
    within.cc <- diag(cc.mean)
    diag(cc.mean) <- 0
    max.between.cc <- apply(cc.mean, 1, max)
    cc.sep <- within.cc - max.between.cc
    cc.stats <- rbind(cc.stats, data.frame(sample.type = sample.type1, 
                                           read.set = sub(paste0(sample.type1, "_"), "", set1), 
                                           cluster = names(within.cc),
                                           within.cc, max.between.cc, cc.sep))
  }
}
cc.stats$read.set <- factor(cc.stats$read.set, 
                            levels = c("IE_IE", "IE_E", "E_IE", "E_E"))
cc.stats$sample.type <- factor(cc.stats$sample.type, 
                            levels = c("nuc", "cell"))



#### Figure 3C - Cluster cohesion/separation summary ####
#getting the convex hull of each unique point set
find_hull <- function(df) df[chull(df$within.cc, df$cc.sep), ]
cc.hulls <- ddply(cc.stats, c("sample.type", "read.set"), find_hull)
# cc.hulls <- subset(cc.hulls, read.set %in% c("IE_IE", "E_E"))

g1 <- ggplot(cc.stats, aes(x = within.cc, y = cc.sep, 
                           color = sample.type, fill = sample.type, 
                           shape = sample.type)) +
  facet_wrap(~ read.set, ncol = 2) +
  geom_polygon(data = cc.hulls, alpha = 0.1, color = NA) +
  geom_point(size = 2) +
  xlim(c(0.7, 1)) +
  xlab("Cluster cohesion") +
  ylab("Cluster separation") +
  scale_color_brewer(palette="Set1") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        strip.background = element_blank(),
        strip.text.x = element_blank())
plot(g1)




#### Within/between co-clustering - paired by matching cluster ####
for (set1 in c("nuc_IE_IE", "cell_IE_IE")) {
  samp.dat <- anno.list[[set1]]
  cl.cons.ordered <- cc.list[[set1]]
  diag(cl.cons.ordered) <- NA
  
  max.mean.between.cc.all <- NULL
  for (i in 1:nrow(samp.dat)) {
    cl.cons1 <- as.numeric(cl.cons.ordered[i, ])
    samp.within <- which(samp.dat$cluster_label == samp.dat$cluster_label[i])
    samp.between <- which(samp.dat$cluster_label != samp.dat$cluster_label[i])
    max.mean.between.cc <- max(tapply(cl.cons1[samp.between], 
                                      samp.dat$cluster_label[samp.between], mean))
    max.mean.between.cc.all <- c(max.mean.between.cc.all, max.mean.between.cc)
  }
  anno.list[[set1]]$cc.mean.max.between <- max.mean.between.cc.all
  anno.list[[set1]]$cc.mean.min.diff <- anno.list[[set1]]$cc.mean.within_label - anno.list[[set1]]$cc.mean.max.between
}



keep.anno.cols <- c("cluster_label", "cc.mean.within_label", 
                    "cc.mean.max.between", "cc.mean.min.diff")
anno.nuc.cell <- rbind(data.frame(sample.type = "nuc", anno.list[["nuc_IE_IE"]][, keep.anno.cols]),
                       data.frame(sample.type = "cell", anno.list[["cell_IE_IE"]][, keep.anno.cols]))
colnames(anno.nuc.cell) <- sub("_label", "", colnames(anno.nuc.cell))
anno.nuc.cell$sample.type <- factor(anno.nuc.cell$sample.type, levels = c("nuc", "cell"))
anno.nuc.cell$cluster_pair <- sapply(anno.nuc.cell$cluster, 
                                     function(x) strsplit(as.character(x), "_")[[1]][4])             

# Cell v nuc cocl metrics
cellvnuc.list <- list()
for (cc.metric in c("cc.mean.within", 
                    "cc.mean.max.between", 
                    "cc.mean.min.diff")) {
  cellvnuc.cl <- NULL
  for (pair1 in unique(anno.nuc.cell$cluster_pair)) {
    anno.nuc.cell.subset <- subset(anno.nuc.cell, cluster_pair == pair1)
    lm1 <- lm(as.formula(paste(cc.metric, "~ sample.type")), anno.nuc.cell.subset)
    coef1 <- t(coef(summary(lm1))["sample.typenuc", ])
    cellvnuc.cl <- rbind(cellvnuc.cl, 
                         data.frame(cluster_pair = pair1, coef1))
  }
  cellvnuc.cl$p_adj <- p.adjust(cellvnuc.cl$Pr...t.., "bonferroni")
  cellvnuc.list[[cc.metric]] <- cellvnuc.cl
}

# lm2 <- lm(cc.mean.min.diff ~ sample.type + cluster_pair, anno.nuc.cell)
# summary(lm2)


# Cell v nuc cocl plots
g2a <- ggplot(anno.nuc.cell, aes(x = cluster_pair, y = cc.mean.within, 
                                 color = sample.type)) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15), 
             size = 0.5, alpha = 0.1) +
  geom_boxplot(outlier.shape = NA, size = 0.3) +
  xlab("Matched clusters") +
  ylab("Cluster cohesion") +
  scale_color_brewer(palette="Set1", name = "", labels = c("Nuclei", "Cells")) +
  theme_bw() +
  theme(panel.grid.major = element_line(size = 0.1), 
        panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45,hjust=1))
plot(g2a)
ggsave(g2a, width = 6, height = 4, 
       filename = "output/nuc_vs_cell_cc.mean.within.pdf")


g2b <- ggplot(anno.nuc.cell, aes(x = cluster_pair, y = cc.mean.max.between, 
                                 color = sample.type)) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15), 
             size = 0.5, alpha = 0.1) +
  geom_boxplot(outlier.shape = NA, size = 0.3) +
  xlab("Matched clusters") +
  ylab("Cluster relatedness") +
  scale_color_brewer(palette="Set1", name = "", labels = c("Nuclei", "Cells")) +
  theme_bw() +
  theme(panel.grid.major = element_line(size = 0.1), 
        panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45,hjust=1))
plot(g2b)
ggsave(g2b, width = 6, height = 4, 
       filename = "output/nuc_vs_cell_cc.mean.max.between.pdf")


#### Figure 4E - Nuclei vs cell cluster separation ####
g2c <- ggplot(anno.nuc.cell, aes(x = cluster_pair, y = cc.mean.min.diff, 
                                 color = sample.type)) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15), 
             size = 0.5, alpha = 0.1) +
  geom_boxplot(outlier.shape = NA, size = 0.3) +
  xlab("Matched clusters") +
  ylab("Cluster separation") +
  scale_color_brewer(palette="Set1", name = "", labels = c("Nuclei", "Cells")) +
  theme_bw() +
  theme(panel.grid.major = element_line(size = 0.1), 
        panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45,hjust=1))
plot(g2c)
ggsave(g2c, width = 6, height = 4, 
       filename = "output/nuc_vs_cell_cc.mean.min.diff.pdf")


