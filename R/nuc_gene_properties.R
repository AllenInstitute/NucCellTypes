# Load libraries
library(plyr)

#### Calc nuclear fraction by type - Top 3 genes ####
load("//allen/programs/celltypes/workgroups/hct/CT_clustering/mouse_cell_vs_nuc/scriptsAndOutputJAM/outputFiles_new/info.RData", verbose = TRUE)

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
anno.nuc <- as.data.frame(read_feather("//allen/programs/celltypes/workgroups/hct/CT_clustering/mouse_cell_vs_nuc/shinydb/20170818_VISp_L5_nuc_1/anno.feather"))
anno.cell <- as.data.frame(read_feather("//allen/programs/celltypes/workgroups/hct/CT_clustering/mouse_cell_vs_nuc/shinydb/20170818_VISp_L5_cell_1/anno.feather"))


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
write.csv(nuc.ratio.df, file = "analysis/nuc_vs_cell/nuc.ratio.estimate.csv",
          row.names = FALSE)


#### Plot comparison of nuclear fraction estimates ####
g.nucfrac <- ggplot(nuc.ratio.df, aes(x = expr_ratio, y = intron_ratio, label = cell_type)) +
  # geom_abline(intercept = 0, slope = 1) +
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
ggsave(g.nucfrac, filename = "analysis/nuc_vs_cell/nuc_frac_by_type.pdf", 
       width = 3, height = 3)

cor(nuc.ratio.df$intron_ratio, nuc.ratio.df$expr_ratio)
lm.2p <- lm(intron_ratio ~ 1 + expr_ratio, nuc.ratio.df)
lm.1p <- lm(intron_ratio ~ 0 + expr_ratio, nuc.ratio.df)
anova(lm.2p, lm.1p)
coef(summary(lm.1p))


#### Plot top 3 gene expression ####
# Load anno
shared.anno.cols <- intersect(colnames(anno.nuc), colnames(anno.cell))
anno.all <- rbind(anno.nuc[, shared.anno.cols], anno.cell[, shared.anno.cols])
anno.all$cell_prep_type_label <- factor(anno.all$cell_prep_type_label, 
                                        levels = c("Nuclei", "Cells"))
anno.all$cell_type <- sapply(anno.all$cluster_label, 
                             function(x) strsplit(x, "_")[[1]][4])
anno.all$cell_type <- factor(anno.all$cell_type, levels = rev(nuc.ratio.df$cell_type))

# Load expression data
expr.nuc <- as.data.frame(read_feather("//allen/programs/celltypes/workgroups/hct/CT_clustering/mouse_cell_vs_nuc/shinydb/20170818_VISp_L5_nuc_exon_1/data_t.feather"))
row.names(expr.nuc) <- expr.nuc$gene
expr.nuc <- expr.nuc[, -1]

expr.cell <- as.data.frame(read_feather("//allen/programs/celltypes/workgroups/hct/CT_clustering/mouse_cell_vs_nuc/shinydb/20170818_VISp_L5_cell_exon_1/data_t.feather"))
row.names(expr.cell) <- expr.cell$gene
expr.cell <- expr.cell[, -1]

nuc.plot.genes <- c("Malat1", "Meg3", "Snhg11")
expr.anno.df <- data.frame(rbind(t(expr.nuc[nuc.plot.genes, ]), 
                                 t(expr.cell[nuc.plot.genes, ])), 
                           anno.all[, c("cell_type", "cell_prep_type_label")])


for (gene1 in nuc.plot.genes) {
  g.expr.box <- ggplot(expr.anno.df, aes(x = cell_type, 
                                         y = log2(get(gene1) + 1),
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
  ggsave(g.expr.box, filename = paste0("analysis/nuc_vs_cell/", 
                                       gene1, "_expr_by_type.pdf"), 
         width = 4, height = 3)
}



#### Plot intronic read comparison ####
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
ggsave(g.intron.box, filename = "analysis/nuc_vs_cell/intron_reads_by_type.pdf", 
       width = 4, height = 3)











#### Nuc fraction individual genes - properties ####
nuc.frac <- read.csv("C:/Users/trygveb/Dropbox/AIBS/Transcriptomics/Manuscripts/WholeCell_vs_Nuc/Tables/nuclearFraction_byCellType.csv", row.names = 1)
nuc.frac$ScaledAverage[which(apply(nuc.frac[, 1:11], 1, function(x) all(is.na(x))))] <- NA
colnames(nuc.frac) <- sub("ScaledAverage", "nucfrac", colnames(nuc.frac))

nuc.marker <- read.csv("C:/Users/trygveb/Dropbox/AIBS/Transcriptomics/Manuscripts/WholeCell_vs_Nuc/Tables/nuc.marker.scores.csv")
cell.marker <- read.csv("C:/Users/trygveb/Dropbox/AIBS/Transcriptomics/Manuscripts/WholeCell_vs_Nuc/Tables/cell.marker.scores.csv")
gene.info <- read.csv("//allen/programs/celltypes/workgroups/humancelltypes/Analyses/CellTypes/Cell_vs_Nuc/Cell_vs_Nuc_manuscript/lib/Mus_musculus.gene_info.csv", na.strings = "-")

all.gene.info <- merge(nuc.marker, cell.marker, by = "gene", suffixes = c("_nuc", "_cell"))
all.gene.info <- merge(all.gene.info, nuc.frac, by.x = "gene", by.y = "row.names")

# Match gene info
# gene.syn <- lapply(gene.info$Synonyms, function(x) strsplit(as.character(x), "|", fixed = "TRUE")[[1]])
match.gene.ids <- match(all.gene.info$gene, gene.info$Symbol)
match.entrez <- match(sub("^LOC", "", all.gene.info$gene[is.na(match.gene.ids)]), 
                      gene.info$GeneID)
match.gene.ids[is.na(match.gene.ids)] <- match.entrez
match.syn <- sapply(all.gene.info$gene[is.na(match.gene.ids)], 
               function(x) {
                 gene.idx <- which(grepl(paste0("^", x, "$"), gene.info$Synonyms) |
                                     grepl(paste0("^", x, "\\|"), gene.info$Synonyms) |
                                     grepl(paste0("\\|", x, "\\|"), gene.info$Synonyms) |
                                     grepl(paste0("\\|", x, "$"), gene.info$Synonyms))
                 unique.idx <- ifelse(length(gene.idx) == 1, gene.idx, NA)
               })
match.gene.ids[is.na(match.gene.ids)] <- match.syn
all.gene.info <- data.frame(all.gene.info, gene.info[match.gene.ids, ])

# Save gene info
save.cols <- c("gene", "GeneID", "Synonyms", "chromosome", "description", "type_of_gene",
               "inh.clusters_nuc", "exc.clusters_nuc",
               "fpkm.max.cluster_nuc", "fpkm.max_nuc", "beta.prop1_nuc",
               "inh.clusters_cell", "exc.clusters_cell",
               "fpkm.max.cluster_cell", "fpkm.max_cell", "beta.prop1_cell", 
               "nucfrac")
all.gene.info.tosave <- all.gene.info[, save.cols]
colnames(all.gene.info.tosave) <- sub("beta.prop1", "marker.score", 
                                      colnames(all.gene.info.tosave), fixed = TRUE)
colnames(all.gene.info.tosave) <- sub("nucfrac", "nuclear.prop", 
                                      colnames(all.gene.info.tosave), fixed = TRUE)
write.csv(all.gene.info.tosave, file = "C:/Users/trygveb/Dropbox/AIBS/Transcriptomics/Manuscripts/WholeCell_vs_Nuc/Tables/gene_info.csv", row.names = FALSE)


# Process gene info
all.gene.info$nucfrac.bin <- cut(all.gene.info$nucfrac, seq(0, 1, 0.1), include.lowest = TRUE)
levels(all.gene.info$nucfrac.bin) <- paste0(seq(0, 0.9, 0.1), "-", seq(0.1, 1, 0.1))
all.gene.info.subset <- subset(all.gene.info, nucfrac <= 1 &
                                 ((total.clusters_nuc > 0 & fpkm.max_nuc > 1) | 
                                    (total.clusters_cell > 0 & fpkm.max_cell > 1)))
all.gene.info.subset$type_of_gene <- mapvalues(all.gene.info.subset$type_of_gene, 
                                               from = c("ncRNA", "protein-coding", "pseudo"),
                                               to = c("Non-coding", "Protein-coding", "Pseudogene"))



# Nuclear enriched cell type markers
paste(subset(all.gene.info.subset, beta.prop1_nuc > 0.4 & nucfrac > 0.8 & fpkm.max_nuc > 32)$gene, collapse = ",")

# Cytoplasm enriched cell type markers
paste(subset(all.gene.info.subset, beta.prop1_cell > 0.4 & nucfrac < 0.05 & fpkm.max_cell > 32)$gene, collapse = ",")


gene.types <- c("Non-coding", "Protein-coding", "Pseudogene")
g.nucfrac.hist <- ggplot(subset(all.gene.info.subset, type_of_gene %in% gene.types), 
                         aes(x = nucfrac)) +
  facet_wrap(~ type_of_gene, scale = "free_y") +
  geom_histogram(color = "black", fill = "grey", lwd = 0.1) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  xlab("Estimated nuclear proportion of transcripts") +
  ylab("Number of genes") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12), 
        panel.border = element_rect(colour = "black"))
plot(g.nucfrac.hist)
ggsave(g.nucfrac.hist, filename = "analysis/nuc_vs_cell/nucfrac_hist_by_genetype.pdf", 
       width = 6, height = 2)



# Cytoplasm enriched genes are less specifically expressed (majority are in all clusters)
g.marker.box <- ggplot(subset(all.gene.info.subset, type_of_gene %in% gene.types), 
                       aes(x = nucfrac.bin, y = beta.prop1_cell)) +
  facet_wrap(~ type_of_gene) +
  geom_boxplot(fill = "grey", outlier.color = "grey90",
               outlier.size = 0.2, lwd = 0.1, fatten = 8) +
  xlab("Estimated nuclear proportion of transcripts") +
  ylab("Marker score") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12), 
        panel.border = element_rect(colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot(g.marker.box)
ggsave(g.marker.box, filename = "analysis/nuc_vs_cell/marker_score_boxplot_by_genetype.pdf", 
       width = 6, height = 3)

# ANOVA
for (type1 in gene.types) {
  aov1 <- aov(beta.prop1_cell ~ nucfrac.bin, 
              data = subset(all.gene.info.subset, type_of_gene == type1))
  pw.diff1 <- TukeyHSD(aov1)
  print(summary(aov1))
  print(pw.diff1[[1]][pw.diff1[[1]][, 4] < 0.05 / 3, ])
}


# Save pseudogenes
write.csv(subset(all.gene.info.subset, nucfrac < 0.05 & type_of_gene == "Pseudogene")$gene, 
          "analysis/nuc_vs_cell/psegenes.csv")  # Cyto
write.csv(subset(all.gene.info.subset, nucfrac > 0.5 & type_of_gene == "Pseudogene")$gene, 
          "analysis/genes.csv")  # Nuclear



#### Compare nuclear proportion to literature (Halpern et al. 2015) ####
nucprop.halpern <- read.csv("//allen/programs/celltypes/workgroups/humancelltypes/Analyses/CellTypes/Cell_vs_Nuc/Cell_vs_Nuc_manuscript/lib/Halpern2015-nuc_cyto_localization/Halpern2015_TableS2_Nuc_Cyto_gene_counts.csv")

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


# Plot comparison of nuc proportions - scatter
cor1 <- round(cor(all.gene.info.subset2$nucfrac, 
                  all.gene.info.subset2$nucprop_mean, use = "pair"), 2)
lm1 <- lm(all.gene.info.subset2$nucfrac ~ 0 + all.gene.info.subset2$nucprop_mean)
slope1 <- round(coef(summary(lm1))[1], 2)

g.nucprop.compare <- ggplot(all.gene.info.subset2, 
                            aes(x = nucfrac, y = nucprop_mean)) +
  geom_point(alpha = 0.1) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method='lm', se = FALSE, formula = y ~ 0 + x, 
              color = "light blue", fullrange = TRUE, size = 1) +
  geom_point(data = subset(all.gene.info.subset2, gene %in% c("Malat1", "Meg3", "Snhg11")),
             aes(x = nucfrac, y = nucprop_mean)) +
  geom_text_repel(data = subset(all.gene.info.subset2, gene %in% c("Malat1", "Meg3", "Snhg11")),
                  aes(x = nucfrac, y = nucprop_mean, label = gene, fontface = "italic")) +
  xlab("Estimated nuclear proportion (this study)") +
  ylab("Estimated nuclear proportion (Halpern et al. 2015)") +
  ggtitle(paste0("Correlation = ", cor1, "; Slope = ", slope1)) +
  coord_fixed() +
  theme_bw()
plot(g.nucprop.compare)
ggsave(g.nucprop.compare, width = 4, height = 4,
       filename = "analysis/nuc_vs_cell/nucprop_compare_halpern2015.pdf")


# Plot comparison of nuc proportions - hist
distrib.diff <- ks.test(all.gene.info.subset2$nucfrac, all.gene.info.subset2$nucprop_mean)

all.gene.info.subset2l <- melt(all.gene.info.subset2[, c("gene", "nucfrac", "nucprop_mean")],
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
       filename = "analysis/nuc_vs_cell/nucprop_compare_hist_halpern2015.pdf")






#### Extra code ####

g1 <- ggplot(data = subset(all.gene.info.subset, beta.prop1_cell > 0 & nucfrac < 1), 
             aes(x = beta.prop1_cell, y = nucfrac)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  theme_bw()
plot(g1)


pheatmap(cor(all.gene.info.subset[, c(2:4, 7:14, 17:24)], use = "pair"))





boxplot(nucfrac ~ type_of_gene, 
        data = all.gene.info.subset,varwidth = TRUE, 
        ylab = "Nuclear fraction",
        col = "grey", las = 2)

boxplot(beta.prop1_cell ~ type_of_gene,
        data = all.gene.info.subset,varwidth = TRUE, 
        ylab = "Marker score (cell)",
        ylim = c(0, 1), las = 2)
boxplot(beta.prop1_nuc ~ type_of_gene, 
        data = all.gene.info.subset,varwidth = TRUE, 
        ylab = "Marker score (nuc)",
        ylim = c(0, 1), las = 2)



par(mfrow = c(1, 3))
for (type1 in gene.types) {
  aov1 <- aov(beta.prop1_cell ~ nucfrac.bin, 
              data = subset(all.gene.info.subset,  type_of_gene == type1))
  print(summary(aov1))
  print(TukeyHSD(aov1))
  
  boxplot(beta.prop1_cell ~ nucfrac.bin, 
          data = subset(all.gene.info.subset,  type_of_gene == type1), 
          varwidth = TRUE, col = "grey", las = 2, main = paste("cell -", type1), 
          xlab = "Nuclear fraction", ylab = "Marker score")
}


aov1 <- aov(beta.prop1_cell ~ type_of_gene, 
            data = subset(all.gene.info.subset, type_of_gene %in% gene.types))
summary(aov1)
TukeyHSD(aov1)

g1b <- ggplot(subset(all.gene.info.subset, type_of_gene %in% gene.types), 
              aes(x = type_of_gene, y = beta.prop1_cell)) +
  geom_boxplot(fill = "grey", outlier.color = "grey90") +
  xlab("Class of gene") +
  ylab("Marker score") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot(g1b)


g1c <- ggplot(all.gene.info, 
              aes(x = type_of_gene, y = nucfrac)) +
  geom_boxplot(fill = "grey", outlier.color = "grey90") +
  xlab("Class of gene") +
  ylab("Nuclear prop.") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot(g1c)





#### Estimate overall nuclear fraction ####
type.nuc.frac <- read.csv("C:/Users/trygveb/Dropbox/AIBS/Transcriptomics/Manuscripts/WholeCell_vs_Nuc/Tables/nuclearFraction_snoRNA.csv")
type.nuc.frac$cell_type <- factor(type.nuc.frac$cell_type, 
                                  levels = type.nuc.frac$cell_type[order(type.nuc.frac$nucfrac_mean)])




#### Plot nuc vs. cell expression with nuclear fraction estimate ####
shared.genes <- intersect(row.names(expr.nuc), row.names(expr.cell))
expr.all <- data.frame(expr.nuc[shared.genes, ], expr.cell[shared.genes, ])
expr.mean.nuc <- t(apply(expr.nuc, 1, function(x) tapply(x, anno.nuc$cluster_label, mean)))
expr.mean.cell <- t(apply(expr.cell, 1, function(x) tapply(x, anno.cell$cluster_label, mean)))
expr.mean.all <- data.frame(expr.mean.nuc[shared.genes, ], expr.mean.cell[shared.genes, ])

malat1.med <- log2(tapply(as.numeric(expr.all["Malat1", ]), cl.lab2, median) + 1)
malat1.diff <- malat1.med[seq(1, length(malat1.med), 2)] - malat1.med[seq(2, length(malat1.med), 2)]


# Plot single type
type1 <- "Samd3"
plot.types <- colnames(expr.mean.all)[grep(type1, colnames(expr.mean.all))]
nuc.frac1 <- type.nuc.frac$nucfrac_mean[grep(type1, type.nuc.frac$cell_type)]
genes.expr1 <- shared.genes[apply(expr.mean.all[, plot.types], 1, function(x) all(x > 0.1))]
nuc.loc.genes <- c(gene.info$Symbol[gene.info$type_of_gene == "snoRNA"],
                   "Malat1", "Meg3", "Snhg11")
plot.nuc.genes <- intersect(genes.expr1, nuc.loc.genes)
plot(log2(expr.mean.all[plot.nuc.genes, plot.types]))
text(log2(expr.mean.all[plot.nuc.genes, plot.types]), labels = plot.nuc.genes, cex = 0.7)
abline(log2(nuc.frac1), 1)

plot(log2(expr.mean.all[, plot.types] + 1), main = type1)
points(log2(expr.mean.all[plot.nuc.genes, plot.types]),
       col = "red", pch = 18, cex = 2)
abline(lm(expr.mean.all[plot.nuc.genes, plot.types[1]] ~ expr.mean.all[plot.nuc.genes, plot.types[2]]))
abline(log2(0.1), 1)
abline(log2(0.3), 1)


g4 <- ggplot(expr.mean.all, aes(x = get(plot.types[1]) + 1, 
                                y = get(plot.types[2]) + 1)) +
  geom_point(alpha = 0.1) +
  geom_abline(intercept = log10(nuc.frac1), slope = 1) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw()
plot(g4)

