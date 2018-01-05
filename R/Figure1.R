library(pheatmap)
library(RColorBrewer)
library(ggplot2)


outpath = "output"

# Load cell/nucleus matching pairs
load(file = "data/map_nuc_to_cells.rda")

#### Figure 1B ####
pdf(file = paste0(out.path, "/map_nuc_to_cells.pdf"), height = 6, width = 6)
plot(max.cor.nuc, max.cor.all, xlab = "Best matching nucleus",
     ylab = "Best matching cell", main = "Expression correlation")
abline(a = 0, b = 1)



#### SI Figure 1A,B ####
# Cre-line / layer summary
cell.subset <- which(samp.cell$exp_component_name %in% samp.ids[["cell"]])
cre.info <- as.character(samp.cell$transgenic_recombinase[cell.subset])
cre.info[cre.info == ""] <- "Virally labeled"
cre.info[which(samp.cell$exp_component_name[cell.subset] == "SM-D9E66_S95_E1-50")] <- "Sim1"
cre.info <- sub("Rorb_neo", "Rorb", cre.info)
cre.info <- sub("Chrna2_Pvalb-T2A-Dre", "Chrna2_Pvalb", cre.info)
cre.info <- sub("Pvalb-IRES-Cre", "Pvalb", cre.info)
cre.cnt <- table(cre.info)

roi.info <- as.character(samp.cell$roi[cell.subset])
roi.info <- sub("VISp_", "", roi.info)
roi.info <- sub("L1-L2-3-L4-L5-L6", "L1-6", roi.info)
roi.info <- sub("L1-L2-3-L4", "L1-4", roi.info)
roi.info <- sub("L2-3-L4", "L2-4", roi.info)
roi.info <- sub("L4-L5-L6", "L4-6", roi.info)
roi.info <- sub("L4-L5", "L4-5", roi.info)
roi.info <- sub("L5-L6", "L5-6", roi.info)
roi.cnt <- table(roi.info)
roi.cnt <- roi.cnt[c(1, 4, 2, 5, 3, 6:11)]


pdf(file = paste0(out.path, "/mapped_cell_cre.pdf"), height = 5, width = 4)
par(mar = c(5, 7, 3, 3))
barplot(sort(cre.cnt / sum(cre.cnt), decreasing = TRUE), 
        xlab = "Proportion of matched cells", xlim = c(0, 0.3),
        main = "Mouse Cre-line",
        col = "grey", horiz = TRUE, las = 1)
dev.off()

pdf(file = paste0(out.path, "/mapped_cell_layer.pdf"), height = 5, width = 4)
par(mar = c(5, 7, 3, 3))
barplot(sort(roi.cnt / sum(roi.cnt), decreasing = TRUE), 
        xlab = "Proportion of matched cells", xlim = c(0, 0.3),
        main = "Dissected layer(s)",
        col = "grey", horiz = TRUE, las = 1)
dev.off()

