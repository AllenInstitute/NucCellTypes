# Load libraries
library(plyr)
library(feather)
library(ggplot2)
library(ggrepel)

#### Calc nuclear fraction by type - Top 3 genes ####
load("../data/expr_summary.rda", verbose = TRUE)

top3.genes <- c("Malat1", "Meg3", "Snhg11")
nuc.ratio.df <- matrix(NA, ncol(meansIC), 6,
                       dimnames = list(colnames(meansIC), 
                                       paste0(top3.genes, c(rep("_mean", 3), rep("_sd", 3)))))
for (gene1 in top3.genes) {
  ratio.mean <- meansIC[gene1, ] / meansIN[gene1, ]
  ratio.var <- ratio.mean^2 * ((sdsIN[gene1, ] / meansIN[gene1, ])^2 +
                                 (sdsIC[gene1, ] / meansIC[gene1, ])^2)
  nuc.ratio.df[, paste0(gene1, "_mean")] <- ratio.mean
  nuc.ratio.df[, paste0(gene1, "_sd")] <- sqrt(ratio.var)
}

ratio.mean <- apply(nuc.ratio.df, 1, function(x) mean(x[1:3]))
ratio.var <- apply(nuc.ratio.df, 1, function(x) sum((x[4:6] / 3)^2))
nuc.ratio.df <- data.frame(cell_type = row.names(nuc.ratio.df), nuc.ratio.df, 
                           expr_ratio = ratio.mean, expr_ratio_sd = sqrt(ratio.var))



#### Calc nuclear fraction by type - Introns ####
# Load annotation data
anno.nuc <- as.data.frame(read_feather("../data/20170818_VISp_L5_nuc/anno.feather"))
anno.cell <- as.data.frame(read_feather("../data/20170818_VISp_L5_cell/anno.feather"))
shared.anno.cols <- intersect(colnames(anno.nuc), colnames(anno.cell))
anno.all <- rbind(anno.nuc[, shared.anno.cols], anno.cell[, shared.anno.cols])
anno.all$cell_prep_type_label <- factor(anno.all$cell_prep_type_label, 
                                        levels = c("Nuclei", "Cells"))
anno.all$cell_type <- sapply(anno.all$cluster_label, 
                             function(x) strsplit(x, "_")[[1]][4])
anno.all$cell_type <- factor(anno.all$cell_type, levels = rev(nuc.ratio.df$cell_type))


nuc.intron <- tapply(anno.nuc$percent_reads_aligned_intron_label, anno.nuc$cluster_label, mean)
nuc.intron.sd <- tapply(anno.nuc$percent_reads_aligned_intron_label, anno.nuc$cluster_label, sd)
cell.intron <- tapply(anno.cell$percent_reads_aligned_intron_label, anno.cell$cluster_label, mean)
cell.intron.sd <- tapply(anno.cell$percent_reads_aligned_intron_label, anno.cell$cluster_label, sd)

nuc.ratio.df <- data.frame(nuc.ratio.df, intron_ratio = NA, intron_ratio_sd = NA)
for (i in 1:nrow(nuc.ratio.df)) {
  nuc.idx <- grep(row.names(nuc.ratio.df)[i], names(nuc.intron))
  cell.idx <- grep(row.names(nuc.ratio.df)[i], names(cell.intron))
  intron.ratio <- cell.intron[cell.idx] / nuc.intron[nuc.idx]
  intron.ratio.var <- intron.ratio^2 * ((cell.intron.sd[cell.idx] / cell.intron[cell.idx])^2 +
                                          (nuc.intron.sd[nuc.idx] / nuc.intron[nuc.idx])^2)
  nuc.ratio.df[i, c("intron_ratio", "intron_ratio_sd")] <- c(intron.ratio, sqrt(intron.ratio.var))
}

order.by.nucfrac <- order(rowMeans(nuc.ratio.df[, c("expr_ratio", "intron_ratio")]))
nuc.ratio.df <- nuc.ratio.df[order.by.nucfrac, ]
write.csv(nuc.ratio.df, file = "../output/nuc.ratio.estimate.csv",
          row.names = FALSE)


# Load expression data
load("../data/20170818_VISp_L5_nuc_exon/20170818_VISp_L5_nuc_exon_iter_cl_data.rda")
expr.nuc <- nbt.data
load("../data/20170818_VISp_L5_cell_exon/20170818_VISp_L5_cell_exon_iter_cl_data.rda")
expr.cell <- nbt.data

nuc.plot.genes <- c("Malat1", "Meg3", "Snhg11")
expr.anno.df <- data.frame(rbind(t(expr.nuc[nuc.plot.genes, ]), 
                                 t(expr.cell[nuc.plot.genes, ])), 
                           anno.all[, c("cell_type", "cell_prep_type_label")])



#### Figure 5A - Plot intronic read comparison ####
g.intron.box <- ggplot(anno.all, aes(x = cell_type, 
                                     y = percent_reads_aligned_intron_label,
                                     color = cell_prep_type_label)) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15), 
             size = 0.5, alpha = 0.1) +
  geom_boxplot(outlier.shape = NA, size = 0.3) +
  scale_color_brewer(palette="Set1", name = "", labels = c("Nuclei", "Cells")) +
  xlab("") +
  ylab("Intron read percentage") +
  ggtitle("Intron read percentage") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot(g.intron.box)
ggsave(g.intron.box, filename = "../output/intron_reads_by_type.pdf", 
       width = 4, height = 3)



#### Figure 5B, S5A - Plot top 3 gene expression ####
for (gene1 in nuc.plot.genes) {
  g.expr.box <- ggplot(expr.anno.df, aes(x = cell_type, y = get(gene1),
                                         color = cell_prep_type_label)) +
    geom_point(position=position_jitterdodge(jitter.width = 0.15), 
               size = 0.5, alpha = 0.1) +
    geom_boxplot(outlier.shape = NA, size = 0.3) +
    scale_color_brewer(palette="Set1", name = "", labels = c("Nuclei", "Cells")) +
    xlab("") +
    ylab("log2(CPM + 1)") +
    ggtitle(paste(gene1, "expression")) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  plot(g.expr.box)
  g.expr.box.fn <- paste0("../output/", gene1, "_expr_by_type.pdf")
  ggsave(g.expr.box, filename = g.expr.box.fn, width = 4, height = 3)
}


#### Figure 5C - Plot comparison of nuclear fraction estimates ####
g.nucfrac <- ggplot(nuc.ratio.df, aes(x = expr_ratio, y = intron_ratio, label = cell_type)) +
  geom_smooth(method='lm', se = FALSE, formula = y ~ 0 + x, 
              fullrange = TRUE, color = "light blue") +
  geom_errorbarh(aes(xmin = expr_ratio - expr_ratio_sd, 
                     xmax = expr_ratio + expr_ratio_sd,
                     height = 0), color = "grey90") +
  geom_errorbar(aes(ymin = intron_ratio - intron_ratio_sd, 
                    ymax = intron_ratio + intron_ratio_sd, 
                    width = 0), color = "grey90") +
  geom_point() +
  geom_text_repel(size = 3) +
  xlim(c(0, 0.7)) +
  ylim(c(0, 0.9)) +
  xlab("Nuclear gene expression ratio") +
  ylab("Intronic read ratio") +
  ggtitle("Nuclear fraction varies among neuron types") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot(g.nucfrac)
ggsave(g.nucfrac, filename = "../output/nuc_frac_by_type.pdf", 
       width = 3, height = 3)

cor(nuc.ratio.df$intron_ratio, nuc.ratio.df$expr_ratio)
lm.2p <- lm(intron_ratio ~ 1 + expr_ratio, nuc.ratio.df)
lm.1p <- lm(intron_ratio ~ 0 + expr_ratio, nuc.ratio.df)
anova(lm.2p, lm.1p)
coef(summary(lm.1p))






#### Nuc fraction individual genes - properties ####
all.gene.info <- read.csv("../data/TableS6_Figure5_gene_info.csv")
all.gene.info$total.clusters_nuc <- all.gene.info$inh.clusters_nuc + all.gene.info$exc.clusters_nuc
all.gene.info$total.clusters_cell <- all.gene.info$inh.clusters_cell + all.gene.info$exc.clusters_cell

# Process gene info
nucprop.bins <- seq(0, 1, 0.1)
all.gene.info$nuclear.prop.bin <- cut(all.gene.info$nuclear.prop, nucprop.bins, include.lowest = TRUE)
levels(all.gene.info$nuclear.prop.bin) <- paste0(nucprop.bins[-length(nucprop.bins)], 
                                            "-", nucprop.bins[-1])
all.gene.info.subset <- subset(all.gene.info, nuclear.prop <= 1 &
                                 ((total.clusters_nuc > 0 & fpkm.max_nuc > 1) | 
                                    (total.clusters_cell > 0 & fpkm.max_cell > 1)))
all.gene.info.subset$type_of_gene <- mapvalues(all.gene.info.subset$type_of_gene, 
                                               from = c("ncRNA", "protein-coding", "pseudo"),
                                               to = c("Non-coding", "Protein-coding", "Pseudogene"))



# Nuclear enriched cell type markers
paste(subset(all.gene.info.subset, marker.score_nuc > 0.4 & nuclear.prop > 0.8 & fpkm.max_nuc > 32)$gene, collapse = ",")

# Cytoplasm enriched cell type markers
paste(subset(all.gene.info.subset, marker.score_cell > 0.4 & nuclear.prop < 0.05 & fpkm.max_cell > 32)$gene, collapse = ",")


#### Figure 5E - Nuclear proportion by gene type ####
gene.types <- c("Non-coding", "Protein-coding", "Pseudogene")
g.nuclear.prop.hist <- ggplot(subset(all.gene.info.subset, type_of_gene %in% gene.types), 
                         aes(x = nuclear.prop, fill = type_of_gene)) +
  facet_wrap(~ type_of_gene, scale = "free_y", ncol = 1) +
  geom_histogram(color = "black", lwd = 0.1, show.legend = FALSE) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Nuclear proportion") +
  ylab("Number of genes") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12), 
        panel.border = element_rect(colour = "black"))
plot(g.nuclear.prop.hist)
ggsave(g.nuclear.prop.hist, filename = "../output/nuclear.prop_hist_by_genetype5.pdf", 
       width = 2.5, height = 5)



#### Figure 5F - Non-coding/pseudogenes better markers ####
# gene <- read.csv("C:/Users/trygveb/Dropbox/AIBS/Transcriptomics/Manuscripts/WholeCell_vs_Nuc/Tables/gene_info.csv")
gene2 <- droplevels(subset(all.gene.info, type_of_gene %in% c("protein-coding", "pseudo", "ncRNA") & 
                             nuclear.prop <= 1 &
                             ((total.clusters_nuc > 0 & fpkm.max_nuc > 1) | 
                                (total.clusters_cell > 0 & fpkm.max_cell > 1))))

kt1 <- kruskal.test(gene2$marker.score_cell ~ gene2$type_of_gene)
pw.wt1 <- pairwise.wilcox.test(gene2$marker.score_cell, gene2$type_of_gene, 
                               paired = FALSE, p.adjust.method = "bonf")


g.marker.bygene <- ggplot(subset(all.gene.info.subset, type_of_gene %in% gene.types), 
                          aes(x = type_of_gene, y = marker.score_cell, fill = type_of_gene)) +
  geom_violin(lwd = 0.1, show.legend = FALSE) +
  geom_boxplot(width = 0.0001, outlier.shape = NA, coef = 0, show.legend = FALSE) +
  stat_summary(fun.y=median, geom="point", size=1, show.legend = FALSE) +
  scale_fill_brewer(palette = "Dark2") +
  xlab("") +
  ylab("Marker score") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12), 
        panel.border = element_rect(colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot(g.marker.bygene)
ggsave(g.marker.bygene, filename = "../output/marker_score_vs_genetype3.pdf", 
       width = 2, height = 2.5)



#### Figure 5G - Cytoplasm enriched genes are less likely marker genes ####
nucprop.bins <- seq(0, 1, 0.1)
gene2$nuclear.prop.bin <- cut(gene2$nuclear.prop, nucprop.bins, include.lowest = TRUE)
levels(gene2$nuclear.prop.bin) <- paste0(nucprop.bins[-length(nucprop.bins)], 
                                         "-", nucprop.bins[-1])
kt2 <- kruskal.test(gene2$marker.score_cell ~ gene2$nuclear.prop.bin)
pw.wt2 <- pairwise.wilcox.test(gene2$marker.score_cell, gene2$nuclear.prop.bin, 
                               paired = FALSE, p.adjust.method = "bonf")
summary(lm(marker.score_cell ~ nuclear.prop, 
           data = subset(all.gene.info.subset, type_of_gene %in% gene.types)))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.214165   0.003018  70.963  < 2e-16 ***
#   nuclear.prop     0.040506   0.007241   5.594 2.27e-08 ***
# TukeyHSD(aov(gene2$marker.score_cell ~ gene2$nuclear.prop.bin))


mean.score <- median(subset(all.gene.info.subset, type_of_gene %in% gene.types)$marker.score_cell)
g.marker.bynucprop <- ggplot(subset(all.gene.info.subset, type_of_gene %in% gene.types), 
                             aes(x = nuclear.prop.bin, y = marker.score_cell)) +
  geom_boxplot(fill = "grey", outlier.color = "grey90",
               outlier.size = 0.2, lwd = 0.1, fatten = 8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, 
              aes(group = 1), col = "light blue") +
  xlab("Nuclear proportion") +
  ylab("Marker score") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12), 
        panel.border = element_rect(colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot(g.marker.bynucprop)
ggsave(g.marker.bynucprop, filename = "../output/marker_score_by_nucprop5.pdf", 
       width = 2.5, height = 2.5)




#### Compare nuclear proportion to literature (Halpern et al. 2015; 10.1016/j.celrep.2015.11.036) ####
nucprop.halpern <- read.csv("../data/Halpern2015_TableS2_Nuc_Cyto_gene_counts.csv")

tissue.samps <- c("MIN6.1", "MIN6.2", "liver.1", "liver.2")
for (tissue1 in tissue.samps) {
  nuc1 <- nucprop.halpern[, paste0("Nuc.", tissue1)]
  cyto1 <- nucprop.halpern[, paste0("Cyto.", tissue1)]
  nucprop1 <- nuc1 / (nuc1 + cyto1)
  nucprop1[which(nuc1 < 1 & cyto1 < 1)] <- NA
  nucprop.halpern <- cbind(nucprop.halpern, nucprop1)
  colnames(nucprop.halpern)[ncol(nucprop.halpern)] <- paste0("nucprop.", tissue1)
}
nucprop.cols <- grep("nucprop", colnames(nucprop.halpern))
cor(nucprop.halpern[, nucprop.cols], use = "pair")
nucprop.halpern$nucprop_mean <- apply(nucprop.halpern[, nucprop.cols], 1, 
                                      mean, na.rm = TRUE)
nucprop.halpern <- na.omit(nucprop.halpern)

all.gene.info.subset2 <- merge(all.gene.info.subset, nucprop.halpern, 
                               by.x = "gene", by.y = "Gene")


#### Figure S5E - Plot comparison of nuc proportions - scatter ####
cor1 <- round(cor(all.gene.info.subset2$nuclear.prop, 
                  all.gene.info.subset2$nucprop_mean, use = "pair"), 2)
lm1 <- lm(all.gene.info.subset2$nuclear.prop ~ 0 + all.gene.info.subset2$nucprop_mean)
slope1 <- round(coef(summary(lm1))[1], 2)

g.nucprop.compare <- ggplot(all.gene.info.subset2, 
                            aes(x = nuclear.prop, y = nucprop_mean)) +
  geom_point(alpha = 0.1) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method='lm', se = FALSE, formula = y ~ 0 + x, 
              color = "light blue", fullrange = TRUE, size = 1) +
  geom_point(data = subset(all.gene.info.subset2, gene %in% c("Malat1", "Meg3", "Snhg11")),
             aes(x = nuclear.prop, y = nucprop_mean)) +
  geom_text_repel(data = subset(all.gene.info.subset2, gene %in% c("Malat1", "Meg3", "Snhg11")),
                  aes(x = nuclear.prop, y = nucprop_mean, label = gene, fontface = "italic")) +
  xlab("Estimated nuclear proportion (this study)") +
  ylab("Estimated nuclear proportion (Halpern et al. 2015)") +
  ggtitle(paste0("Correlation = ", cor1, "; Slope = ", slope1)) +
  coord_fixed() +
  theme_bw()
plot(g.nucprop.compare)
ggsave(g.nucprop.compare, width = 4, height = 4,
       filename = "../output/nucprop_compare_halpern2015.pdf")


#### Figure S5F - Plot comparison of nuc proportions - hist ####
distrib.diff <- ks.test(all.gene.info.subset2$nuclear.prop, all.gene.info.subset2$nucprop_mean)

all.gene.info.subset2l <- melt(all.gene.info.subset2[, c("gene", "nuclear.prop", "nucprop_mean")],
                               id = "gene", value.name = "nucprop")

g.nucprop.hist.compare <- ggplot(all.gene.info.subset2l, 
                                 aes(x = nucprop, fill = variable)) +
  geom_histogram(binwidth = 0.02, alpha = 0.5, position = "identity") +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_fill_brewer(palette="Set1", name = "", 
                    labels = c("This study", "Halpern et al. 2015")) +
  xlab("Estimated nuclear proportion") +
  ylab("Number of genes") +
  theme_bw(base_size = 8) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot(g.nucprop.hist.compare)
ggsave(g.nucprop.hist.compare, width = 5, height = 2,
       filename = "../output/nucprop_compare_hist_halpern2015.pdf")
