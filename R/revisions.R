library(feather)

anno <- as.data.frame(read_feather( "//allen/programs/celltypes/workgroups/hct/CT_clustering/mouse_cell_vs_nuc/shinydb/20170818_VISp_L5_nuc_1/anno_intergenic_modified.feather"))

apply(anno[, c("percent_reads_aligned_exon_label", 
               "percent_reads_aligned_intron_label",
               "percent_reads_intergenic_label")], 1, sum)


load("//allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/facs/R_Object/20180710_RSC-11-161_mouse_star_samp.dat.Rdata", verbose = TRUE)



#### Compare intron/exon gene detection ####

# Load data
gene.info <- read.csv(file = "data/TableS6_Figure5_gene_info.csv")

# Exons + Introns
load("data/20170818_VISp_L5_nuc/20170818_VISp_L5_nuc_iter_cl_data.rda", verbose = TRUE)
samp.dat.nuc <- samp.dat
nbt.data.nuc <- nbt.data
load("data/20170818_VISp_L5_cell/20170818_VISp_L5_cell_iter_cl_data.rda", verbose = TRUE)
samp.dat.cell <- samp.dat
nbt.data.cell <- nbt.data

# Exons only
load("data/20170818_VISp_L5_nuc_exon/20170818_VISp_L5_nuc_exon_iter_cl_data.rda", verbose = TRUE)
samp.dat.nuc.exon <- samp.dat
nbt.data.nuc.exon <- nbt.data
load("data/20170818_VISp_L5_cell_exon/20170818_VISp_L5_cell_exon_iter_cl_data.rda", verbose = TRUE)
samp.dat.cell.exon <- samp.dat
nbt.data.cell.exon <- nbt.data

# Combine data
samp.dat <- rbind(samp.dat.nuc, samp.dat.cell)
samp.dat$Samples <- as.factor(samp.dat$cell_prep_type)

nuc.cell.genes <- intersect(row.names(nbt.data.nuc), row.names(nbt.data.cell))
cpm.dat <- cbind(nbt.data.nuc[nuc.cell.genes, ], nbt.data.cell[nuc.cell.genes, ])
rownames(cpm.dat) <- toupper(rownames(cpm.dat))
cpm.exon.dat <- cbind(nbt.data.nuc.exon[nuc.cell.genes, ], nbt.data.cell.exon[nuc.cell.genes, ])
rownames(cpm.exon.dat) <- toupper(rownames(cpm.exon.dat))


detect.prop <- data.frame(Cells_e = apply(cpm.exon.dat[, samp.dat$Samples == "Cells"], 1, 
                                          function(x) sum(x > 0) / length(x)),
                          Cells_ie = apply(cpm.dat[, samp.dat$Samples == "Cells"], 1, 
                                          function(x) sum(x > 0) / length(x)),
                          Nuclei_e = apply(cpm.exon.dat[, samp.dat$Samples == "Nuclei"], 1, 
                                          function(x) sum(x > 0) / length(x)),
                          Nuclei_ie = apply(cpm.dat[, samp.dat$Samples == "Nuclei"], 1, 
                                          function(x) sum(x > 0) / length(x)),
                          Cells_e_expr = apply(cpm.exon.dat[, samp.dat$Samples == "Cells"], 1, 
                                          mean))

apply(detect.prop, 2, function(x) sum(x > 0))


incr.detection <- which(detect.prop$Nuclei_ie - detect.prop$Nuclei_e > 0.25)

pdf(file = "output/nuc_intron_detect.pdf", width = 4, height = 4)
plot(density(gene.info$marker.score_cell[detect.prop$Cells_e_expr > 0]), 
     las = 1, lwd = 2, main = "", xlab = "Cell marker score (beta)")
lines(density(gene.info$marker.score_cell[incr.detection]), col = "red", lwd = 2)

plot(density(detect.prop$Cells_e_expr[detect.prop$Cells_e_expr > 0]), 
     # xlim = c(0, 12), 
     las = 1, lwd = 2, main = "", xlab = "Cell mean expression (log2CPM + 1)")
lines(density(detect.prop$Cells_e_expr[incr.detection]), col = "red", lwd = 2)
dev.off()


gene.type <- cbind(table(gene.info$type_of_gene[detect.prop$Cells_e_expr > 0]),
                   table(gene.info$type_of_gene[incr.detection]))
sweep(gene.type, 2, colSums(gene.type), "/")
mosaicplot(t(gene.type[c(1, 3, 4), ]), shade = TRUE)


write.csv(row.names(detect.prop)[incr.detection], "output/incr.detection.csv", 
          row.names = FALSE)




#### Dendritic localization ####
d.genes <- scan("data/Cajigas2012_S10.txt", "character")
d.id <- which(gene.info$gene %in% d.genes & 
                !is.na(gene.info$nuclear.prop) & 
                gene.info$fpkm.max_cell > 1)

plot(log2(gene.info$fpkm.max_cell + 1)[d.id], gene.info$nuclear.prop[d.id])


plot(density(log2(gene.info$fpkm.max_cell + 1)[gene.info$fpkm.max_cell > 1]))
lines(density(log2(gene.info$fpkm.max_cell[d.id] + 1)), col = "blue")

plot(density(gene.info$nuclear.prop[gene.info$fpkm.max_cell > 1]), ylim = c(0, 3))
lines(density(gene.info$nuclear.prop[d.id]), col = "blue")

median(gene.info$nuclear.prop[gene.info$fpkm.max_cell > 1])
median(gene.info$nuclear.prop[d.id])
