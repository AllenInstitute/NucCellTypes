library(dplyr)
library(ggplot2)


scores.nuc <- read.csv("C:/Users/trygveb/Dropbox/AIBS/Transcriptomics/Manuscripts/WholeCell_vs_Nuc/Tables/nuc.marker.scores.csv")
scores.cell <- read.csv("C:/Users/trygveb/Dropbox/AIBS/Transcriptomics/Manuscripts/WholeCell_vs_Nuc/Tables/cell.marker.scores.csv")

scores <- merge(scores.nuc, scores.cell, by = "gene", 
                suffixes = c("_nuc", "_cell"))


# Select best shared markers of types
as_tibble(scores) %>%
  filter(total.clusters_nuc == 1 & 
           # total.clusters_cell == 1 &
           !grepl("^Gm", gene) & !grepl("^LOC", gene) &
           !grepl("^[0-9]", gene)) %>%
  group_by(gene) %>%
  mutate(min_score = min(beta.expr_cell, beta.expr_nuc)) %>%
  group_by(fpkm.max.cluster_cell) %>%
  arrange(fpkm.max.cluster_cell) %>%
  top_n(1, min_score) %>%
  select(gene, min_score, beta.prop1_cell, beta.prop1_nuc, 
         beta.expr_cell, beta.expr_nuc)


scores.subset <- subset(scores, total.clusters_nuc > 0 & total.clusters_nuc < 11 &
                          total.clusters_cell > 0 & total.clusters_cell < 11)
summary(lm(beta.prop1_cell ~ 0 + beta.prop1_nuc, data = scores.subset))

g1 <- ggplot(scores.subset, aes(x = beta.prop1_nuc, y = beta.prop1_cell)) + 
  # color = log2(fpkm.max_cell + 1))) +
  # facet_wrap(~ total.clusters_cell) +
  geom_point(alpha = 0.1) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method='lm', se = FALSE, formula = y ~ 1 + x, 
              fullrange = TRUE, color = "light blue") +
  # geom_text(data = subset(scores.subset, beta.prop1_cell > 0.6 & beta.prop1_nuc < 0.2),
  #           aes(label = gene), color = "black", size = 3, hjust = 0, nudge_x = 0.01) +
  xlim(c(0, 1)) +
  ylim(c(0, 1)) +
  xlab("Nuclei marker score") +
  ylab("Cells marker score") +
  scale_color_gradientn(colors = grey.colors(10, 0.9, 0, gamma = 1.4)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
plot(g1)
ggsave(g1, width = 2.5, height = 2.5, 
       filename = "analysis/nuc_vs_cell/nuc_vs_cell_marker_scores.pdf")



genes <- subset(scores.subset, beta.prop1_cell < 0.3 & beta.prop1_nuc > 0.5)
genes <- genes[order(genes$beta.prop1_cell - genes$beta.prop1_nuc, decreasing = TRUE), "gene"]
genes <- genes[! grepl("^[0-9]", genes)]
paste(genes, collapse = ",")



table(scores.subset$total.clusters_nuc, scores.subset$total.clusters_cell)






#### Save output ####
keep.genes <- intersect(scores.nuc$gene, scores.cell$gene)
scores.nuc <- scores.nuc[match(keep.genes, scores.nuc$gene), ]
scores.cell <- scores.cell[match(keep.genes, scores.cell$gene), ]

write.csv(scores.nuc, file = "C:/Users/trygveb/Dropbox/AIBS/Transcriptomics/Manuscripts/WholeCell_vs_Nuc/Tables/nuc.marker.scores.csv", row.names = FALSE)
write.csv(scores.cell, file = "C:/Users/trygveb/Dropbox/AIBS/Transcriptomics/Manuscripts/WholeCell_vs_Nuc/Tables/cell.marker.scores.csv", row.names = FALSE)







#### Extra code ####
scoresl <- melt(scores[, c("gene", "fpkm.max_cell", 
                           "tau.prop1_cell", "tau.prop1_nuc")], 
                id = c("gene", "fpkm.max_cell"))


# scores.subset <- subset(scores, beta.prop1_nuc > 0)
scores.subset <- scores[order(scores$fpkm.max_cell), ]
scores.subset$fpkm.max_nuc_bin <- cut(scores.subset$fpkm.max_nuc, 0:15, include.lowest = TRUE)
scores.subset$fpkm.max_cell_bin <- cut(scores.subset$fpkm.max_cell, 0:15, include.lowest = TRUE)
scores.subset$beta.prop1_nuc_bin <- cut(scores.subset$beta.prop1_nuc, 
                                        seq(0, 1, length.out = 6), include.lowest = TRUE)
scores.subset$beta.prop1_cell_bin <- cut(scores.subset$beta.prop1_cell, 
                                         seq(0, 1, length.out = 6), include.lowest = TRUE)


g1 <- ggplot(scores.subset, aes(x = beta.prop1_nuc_bin, y = log2(fpkm.max_nuc + 1))) +
  geom_boxplot(varwidth = TRUE) +
  theme_bw()
plot(g1)



g1 <- ggplot(scores.subset, aes(x = beta.prop1_nuc, y = log2(fpkm.max_nuc + 1))) +
  geom_point() +
  geom_text(data = subset(scores.subset, beta.prop1_nuc > 0.82 | 
                            (beta.prop1_nuc > 0.3 & log2(fpkm.max_nuc) > 9)), 
            aes(label = gene), hjust = 0, nudge_x = 0.01) +
  xlim(c(0, 1)) +
  theme_bw()
plot(g1)


g1 <- ggplot(scores.subset, aes(x = beta.prop1_cell, y = log2(fpkm.max_cell + 1))) +
  geom_point() +
  geom_text(data = subset(scores.subset, beta.prop1_cell > 0.95 | 
                            (beta.prop1_cell > 0.3 & log2(fpkm.max_cell) > 10)), 
            aes(label = gene), hjust = 0, nudge_x = 0.01) +
  theme_bw()
plot(g1)




g3 <- ggplot(scores.subset, aes(y = beta.prop1_cell - beta.prop1_nuc, 
                                # x = log2(fpkm.max_cell + 1))) +
                                x = total.clusters_cell - total.clusters_nuc)) +
  geom_point() +
  # geom_text(data = subset(scores.subset, abs(beta.prop1_nuc - beta.prop1_cell) > 0.5), 
  #           aes(label = gene), color = "black", hjust = 0, nudge_x = 0.01) +
  scale_color_gradientn(colors = grey.colors(10, 0.9, 0, gamma = 1)) +
  theme_bw()
plot(g3)



scores.subset <- subset(scores, (total.clusters_nuc == 1 & tau.prop1_nuc > 0.9) |
                          (total.clusters_cell == 1 & tau.prop1_cell > 0.9))

cor(scores.subset$tau.prop1_nuc, scores.subset$tau.prop1_cell)
plot(scores.subset$tau.prop1_nuc, scores.subset$tau.prop1_cell)
abline(a = 0, b = 1)
abline(a = 0.4, b = 1)
abline(a = -0.4, b = 1)
subset(scores, tau.prop1_nuc > 0.6 & tau.prop1_cell < 0.4)
subset(scores, tau.prop1_nuc < 0.2 & tau.prop1_cell > 0.8)


