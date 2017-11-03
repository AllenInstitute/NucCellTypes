library(pheatmap)
library(RColorBrewer)
library(ggplot2)


input.path <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/STARforLIMS/Mouse/Final_Output/VISp_L5"


#### V1 L5 whole cells ####
for (file1 in c("samp.dat", "introns", "exon", "zero.wt")) { 
  load(paste0(input.path, "/Cell/", file1, "_Mouse_Star_VISp_L5_Cells.Rdata"))
}

keep.cells <- with(samp.dat, 
                   which(roi == "VISp_L5" & 
                           Injection_type == 0 &
                           cell_prep_type == "Cells" & 
                           Exclude.1 == "No"))
samp.cell <- droplevels(samp.dat[keep.cells, ])
expr.cell <- log2(countsR[, keep.cells] + introns[, keep.cells] + 1)
wt.cell <- zero.wt$matw


#### V1 L5 nuclei ####
for (file1 in c("samp.dat", "introns", "exon", "zero.wt")) { 
  load(paste0(input.path, "/Nuc/", file1, "_Mouse_Star_VISp_L5_Nuc.Rdata"))
}

excl.samp <- with(samp.dat, 
                  reads_aligned_to_mrna + reads_aligned_to_genome_only < 500000 |
                    percent_reads_aligned_total < 75 |
                    percent_unique_reads < 50 )
excl.samp[is.na(excl.samp)] <- FALSE
keep.nuc <- with(samp.dat, which(cell_prep_type == "Nuclei" & ! excl.samp))
samp.nuc <- droplevels(samp.dat[keep.nuc, ])
expr.nuc <- log2(countsR[, keep.nuc] + introns[, keep.nuc] + 1)
wt.nuc <- zero.wt$matw

#### Find nuc/cell correlations ####
# cor.all <- matrix(NA, ncol(expr.nuc), ncol(expr.cell), 
#                   dimnames = list(colnames(expr.nuc), colnames(expr.cell)))
#                   
# for (nuc1 in 1:ncol(expr.nuc)) {
#   cor.all[nuc1, ] <- sapply(1:ncol(expr.cell), 
#                             function(i) cor(expr.cell[, i], expr.nuc[, nuc1]))
# }

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
order1 <- order(apply(cor.all, 1, max), decreasing = TRUE)
cor.all.ordered <- cor.all[order1, ]
max.cell.all <- NULL
max.cor.all <- NULL
for (i in 1:nrow(cor.all.ordered)) {
  max.cor <- max(cor.all.ordered[i, ])
  max.cell <- which.max(cor.all.ordered[i, ])
  max.cor.all <- c(max.cor.all, max.cor)
  max.cell.all <- c(max.cell.all, max.cell)
  cor.all.ordered[, max.cell] <- 0
}


#### Plot results
pdf(file = "analysis/20170818/map_nuc_to_cells.pdf", height = 4, width = 8)
hist(max.cor.all, breaks = 20, xlim = c(0, 1), cex.lab = 1.2, 
     xlab = "Expression correlation")
# Best matching Cre-lines
pie(with(droplevels(subset(samp.cell, exp_component_name %in% cell.ids)), 
         table(transgenic_recombinase)), col = brewer.pal(12, "Set3"), cex = 0.8)

g1 <- ggplot(data.frame(x = expr.cell[, max.cell.all[300]], 
                        y = expr.nuc[, order1][, 300]), aes(x=x, y=y)) +
  geom_hex() +
  scale_fill_continuous(low = "gray85", high = "black", trans = "log10") +
  xlab("Whole cell expression") +
  ylab("Nucleus expression") +
  theme_bw()
plot(g1)

cor1 <- cor(cbind(expr.cell[, max.cell.all[c(1, 100, 200, 300, 399)]], 
                  expr.nuc[, order1][, c(1, 100, 200, 300, 399)]))
pheatmap(cor1[1:5, 6:10], cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

samp.ids <- list()
samp.ids[["cell"]] <- samp.cell$exp_component_name[max.cell.all]
samp.ids[["nuc"]] <- samp.nuc$exp_component_name[order1]
# Select random set of Snap25+ cells to compare to nuclei
set.seed(123)
samp.ids[["snap25"]] <- sample(samp.cell$exp_component_name[samp.cell$cre_line == "Snap25"],
                     nnuc, replace = FALSE)

for (id.set in c("cell", "nuc", "snap25")) {
  write.table(samp.ids[[id.set]], file = paste0("analysis/20170818/", id.set, "_ids.txt"),
              row.names = FALSE, col.names = FALSE)
}
save(cor.all.ordered, max.cell.all, max.cor.all, samp.ids, 
     file = "analysis/20170818/map_nuc_to_cells.rda")
# setdiff(cell.ids, samp.dat$exp_component_name)


