library(pheatmap)
library(RColorBrewer)
library(ggplot2)


input.path.cell <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/STARforLIMS/Mouse/Final_Output/Mouse_Final_21913"
input.path.cell.nuc <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/STARforLIMS/Mouse/Final_Output/VISp_L5"
out.path <- "//allen/programs/celltypes/workgroups/humancelltypes/Analyses/CellTypes/Cell_vs_Nuc/Cell_vs_Nuc_manuscript/analysis/match_nuc_cell"


#### V1 L5 whole cells ####
for (file1 in c("samp.dat", "intron", "exon", "zero.wt")) { 
  load(paste0(input.path.cell, "/", file1, "_mouse_star_21913.Rdata"))
}

keep.cell <- which(grepl("VISp", samp.dat$roi))

samp.cell <- droplevels(samp.dat[keep.cell, ])
cd.cell <- introns[, keep.cell] + countsR[, keep.cell]
cpm.cell <- round(sweep(cd.cell, 2, colSums(cd.cell), "/") * 1e6, 0)
expr.cell <- log2(cpm.cell + 1)
wt.cell <- zero.wt$matw


#### V1 L5 nuclei ####
for (file1 in c("samp.dat", "introns", "exon", "zero.wt")) { 
  load(paste0(input.path.cell.nuc, "/Nuc/", file1, "_Mouse_Star_VISp_L5_Nuc.Rdata"))
}

excl.samp <- with(samp.dat, 
                  reads_aligned_to_mrna + reads_aligned_to_genome_only < 500000 |
                    percent_reads_aligned_total < 75 |
                    percent_unique_reads < 50 )
excl.samp[is.na(excl.samp)] <- FALSE
keep.nuc <- with(samp.dat, which(cell_prep_type == "Nuclei" & ! excl.samp))
samp.nuc <- droplevels(samp.dat[keep.nuc, ])
cd.nuc <- introns[, keep.nuc] + countsR[, keep.nuc]
cpm.nuc <- round(sweep(cd.nuc, 2, colSums(cd.nuc), "/") * 1e6, 0)
expr.nuc <- log2(cpm.nuc + 1)
wt.nuc <- zero.wt$matw



#### Find nuc/cell correlations ####
expr.nuc.cell <- cbind(expr.nuc, expr.cell)
nnuc <- ncol(expr.nuc)
ncell <- ncol(expr.cell)
cor.all <- matrix(NA, nnuc, ncell,
                  dimnames = list(colnames(expr.nuc), colnames(expr.cell)))
for (nuc1 in 1:nnuc) {
  print(nuc1)
  for (cell1 in 1:ncell) {
    wt1 <- wt.nuc[, nuc1] * wt.cell[, cell1]
    cor1 <- cov.wt(expr.nuc.cell[, c(nuc1, nnuc + cell1)], wt = wt1, 
                   center = FALSE, cor = TRUE)$cor
    cor.all[nuc1, cell1] <- cor1[1, 2]
  }
}

# Select best matching cells
order.cell <- order(apply(cor.all, 1, max), decreasing = TRUE)
cor.all.ordered <- cor.all[order.cell, ]
max.cell.all <- NULL
max.cor.all <- NULL
for (i in 1:nrow(cor.all.ordered)) {
  max.cor <- max(cor.all.ordered[i, ])
  max.cell <- which.max(cor.all.ordered[i, ])
  max.cor.all <- c(max.cor.all, max.cor)
  max.cell.all <- c(max.cell.all, max.cell)
  cor.all.ordered[, max.cell] <- 0
}


# Select best matching Snap25 cells
cor.snap25 <- cor.all[, which(samp.cell$transgenic_recombinase == "Snap25")]
order.snap25 <- order(apply(cor.snap25, 1, max), decreasing = TRUE)
cor.snap25.ordered <- cor.snap25[order.snap25, ]
max.snap25.all <- NULL
max.cor.snap25 <- NULL
for (i in 1:nrow(cor.snap25.ordered)) {
  max.cor <- max(cor.snap25.ordered[i, ])
  max.snap25 <- which.max(cor.snap25.ordered[i, ])
  max.cor.snap25 <- c(max.cor.snap25, max.cor)
  max.snap25.all <- c(max.snap25.all, max.snap25)
  cor.snap25.ordered[, max.snap25] <- 0
}
match.cell.snap25.order <- match(row.names(cor.all.ordered),
                                 row.names(cor.snap25.ordered))
cor.snap25.ordered <- cor.snap25.ordered[match.cell.snap25.order, ]
max.snap25.all <- max.snap25.all[match.cell.snap25.order]
max.cor.snap25 <- max.cor.snap25[match.cell.snap25.order]


#### Find nuc/nuc correlations ####
nnuc <- ncol(expr.nuc)
cor.nuc <- matrix(NA, nnuc, nnuc,
                  dimnames = list(colnames(expr.nuc), colnames(expr.nuc)))
for (nuc1 in 1:nnuc) {
  print(nuc1)
  for (nuc2 in 1:nnuc) {
    wt1 <- wt.nuc[, nuc1]^2
    cor1 <- cov.wt(expr.nuc[, c(nuc1, nuc2)], wt = wt1, 
                   center = FALSE, cor = TRUE)$cor
    cor.nuc[nuc1, nuc2] <- cor1[1, 2]
  }
}
diag(cor.nuc) <- 0


# Select best matching nuclei
cor.nuc.ordered <- cor.nuc[order.cell, order.cell]
max.nuc.all <- apply(cor.nuc.ordered, 1, which.max)
max.cor.nuc <- apply(cor.nuc.ordered, 1, max)



#### Find matching cell/cell correlations ####
expr.cell.subset <- expr.cell[, names(max.cell.all)]
ncell <- ncol(expr.cell.subset)
cor.cell <- matrix(NA, ncell, ncell,
                   dimnames = list(colnames(expr.cell.subset), 
                                   colnames(expr.cell.subset)))
for (cell1 in 1:ncell) {
  print(cell1)
  for (cell2 in 1:ncell) {
    wt1 <- wt.cell[keep.genes, cell1] * wt.cell[keep.genes, cell2]
    cor1 <- cov.wt(expr.cell.subset[keep.genes, c(cell1, cell2)], wt = wt1, 
                   center = FALSE, cor = TRUE)$cor
    cor.cell[cell1, cell2] <- cor1[1, 2]
  }
}
diag(cor.cell) <- 0



# Sample ids
samp.ids <- list()
samp.ids[["nuc"]] <- samp.nuc$exp_component_name[order.cell]
samp.ids[["cell"]] <- names(max.cell.all)
samp.ids[["snap25_matched"]] <- names(max.snap25.all)
samp.ids[["nuc_matched"]] <- samp.nuc$exp_component_name[max.nuc.all][order.cell]
# Select random set of Snap25+ cells to compare to nuclei
set.seed(123)
samp.ids[["snap25_random"]] <- sample(samp.cell$exp_component_name[samp.cell$cre_line == "Snap25"],
                                      nnuc, replace = FALSE)


for (id.set in c("cell", "nuc", "snap25_matched", "snap25_random", "nuc_matched")) {
  write.table(samp.ids[[id.set]], file = paste0(out.path, "/", id.set, "_ids.txt"),
              row.names = FALSE, col.names = FALSE)
}

if (! file.exists(paste0(out.path, "/map_nuc_to_cells.rda"))) {
  save(cor.all.ordered, max.cell.all, max.cor.all, 
       cor.nuc.ordered, max.nuc.all, max.cor.nuc,
       cor.snap25.ordered, max.snap25.all, max.cor.snap25,
       samp.ids, cor.cell,
       file = paste0(out.path, "/map_nuc_to_cells.rda"))
}


#### Plot results
pdf(file = paste0(out.path, "/map_nuc_to_cells.pdf"), height = 6, width = 6)
hist(max.cor.all, breaks = 20, xlim = c(0, 1), cex.lab = 1.2, 
     xlab = "Expression correlation")
plot(max.cor.nuc, max.cor.all, xlab = "Best matching nucleus",
     ylab = "Best matching cell", main = "Expression correlation")
abline(a = 0, b = 1)

plot(max.cor.snap25, max.cor.all, xlab = "Best matching cell (Snap25)",
     ylab = "Best matching cell", main = "Expression correlation")
abline(a = 0, b = 1)

g1 <- ggplot(data.frame(x = expr.cell[, max.cell.all[440]], 
                        y = expr.nuc[, order.cell][, 440]), aes(x=x, y=y)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point() +
  scale_fill_continuous(low = "gray85", high = "black", trans = "log10") +
  xlab("Whole cell expression") +
  ylab("Nucleus expression") +
  theme_bw()
plot(g1)

g2 <- ggplot(data.frame(x = expr.cell[, max.cell.all[231]], 
                        y = expr.nuc[, order.cell][, 231]), aes(x=x, y=y)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point() +
  scale_fill_continuous(low = "gray85", high = "black", trans = "log10") +
  xlab("Whole cell expression") +
  ylab("Nucleus expression") +
  theme_bw()
plot(g2)

g3 <- ggplot(data.frame(x = expr.cell[, max.cell.all[23]], 
                        y = expr.nuc[, order.cell][, 23]), aes(x=x, y=y)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point() +
  scale_fill_continuous(low = "gray85", high = "black", trans = "log10") +
  xlab("Whole cell expression") +
  ylab("Nucleus expression") +
  theme_bw()
plot(g3)

dev.off()



#### Cre-line / layer summary ####
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









#### Extra code ####
barplot(as.matrix(rev(roi.cnt)), col = grey(c(1,2.5,2.5,3,3.5,4,4.5,5,5,5.5,6) / 6), 
        las = 1)
legend("topleft", fill = rev(grey(c(1,2.5,2.5,3,3.5,4,4.5,5,5,5.5,6) / 6)),
       legend = names(roi.cnt), cex = 0.5, bty = "n")

cre.roi.tbl <- table(roi.info, cre.info)
cre.roi.tbl <- cre.roi.tbl[nrow(cre.roi.tbl):1, order(colSums(cre.roi.tbl))]
barplot(cre.roi.tbl, las = 2,
        col = grey(c(1,2.5,2.5,3,3.5,4,4.5,5,5,5.5,6) / 6))