root.path <- "//allen/programs/celltypes/workgroups"
input.path.cell.nuc <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/STARforLIMS/Mouse/Final_Output/Nuclei_Noisemodel/R_objects/"
analysis.path <- "//allen/programs/celltypes/workgroups/humancelltypes/Analyses/CellTypes/Cell_vs_Nuc/analysis/20170818/match_nuc_cell"

# Load functions
source('//allen/programs/celltypes/workgroups/rnaseqanalysis/Script_Repository/Trygve/src/fMatchRowCol.R')
source("//allen/programs/celltypes/workgroups/rnaseqanalysis/Script_Repository/Trygve/src/celltypes/fMaptoV1Tree.R")


# Noise genes
sex.genes <- scan(file = "//allen/programs/celltypes/workgroups/humancelltypes/Analyses/CellTypes/Cell_vs_Nuc/lib/sex_genes_20160118_GRCm38.p3.txt", "character")
mito.genes <- scan(file = "//allen/programs/celltypes/workgroups/humancelltypes/Analyses/CellTypes/Cell_vs_Nuc/lib/mouse_mito_genes.txt", "character")


for (file1 in c("samp.dat", "intron", "exon", "zero.wt")) { 
  load(paste0(input.path.cell.nuc, file1, "_MouseNuclei6960_star.Rdata"))
}

# Remove sex and mito genes
keep.genes <- which(! rownames(intron) %in% c(sex.genes, mito.genes))
introns <- intron[keep.genes, ]
countsR <- exon[keep.genes, ]


# Keep QC passed cells
# samp.ids <- scan(file = paste0(analysis.path, "/", samp.set, "_ids.txt"),
#                  "character")
# qc.samp <- which(samp.dat$exp_component_name %in% samp.ids)
qc.samp <- which(samp.dat$Exclude.1 == "No")

# Keep samples/genes with information in all matrices
samp.dat <- droplevels(samp.dat[qc.samp, ])

# Create expression matrix
cd <- introns[, qc.samp] + countsR[, qc.samp]
cd2 <- cd
cpm <- round(sweep(cd2, 2, colSums(cd), "/") * 1e6, 0)
cpm <- apply(cpm, 2, function(x) {storage.mode(x) <- 'integer'; x})


#### Map to mouse ####
samp.dat$max.leaf <- apply(log2(cpm + 1), 2, 
                           function(x) MaptoV1Tree(expr = x, 
                                                   genes = rownames(cpm),
                                                   root.path = root.path)$leaf.max)
samp.dat$max.leaf.class <- apply(log2(cpm + 1), 2, 
                                 function(x) MaptoV1Tree(expr = x, 
                                                         genes = rownames(cpm),
                                                         root.path = root.path)$leaf.max.class)


# Data list
dat.list <- MatchRowCol(list(cd, cpm, zero.wt$matw), 
                        df.names = c("counts", "cpm", "zero.data"), 
                        row.caps = FALSE)
samp.dat <- droplevels(samp.dat[dat.list$common.col, ])

# Format gene/cell names for automated script
nbt.data <- as.data.frame(data.matrix(log2(dat.list$cpm + 1)))

# Data for DESeq variable gene selection
cnt.data <- as.data.frame(dat.list$counts)
ercc.data <- NULL

# Data for scde
err.data <- zero.wt$o.ifm[dat.list$common.col, ]
zero.data <- as.data.frame(dat.list$zero.data)

# Add roi to sample names
samp.name <- as.character(samp.dat$transgenic_recombinase)
colnames(nbt.data) <- paste0(colnames(nbt.data), "~", samp.name)
colnames(cnt.data) <- paste0(colnames(cnt.data), "~", samp.name)
colnames(zero.data) <- paste0(colnames(zero.data), "~", samp.name)
row.names(err.data) <- paste0(rownames(err.data), "~", samp.name)

# Batch (for scde correction)
batch.id <- NULL

# Add dummy data sets
ercc.data <- NULL

# Save data
dat.path <- "//allen/programs/celltypes/workgroups/hct/CT_clustering/mouse_cell_vs_nuc/data"
out.fn <- paste0(dat.path, "/20171211_VISp_nuc_iter_cl_data.rda")
if (! file.exists(out.fn)) {
  save(nbt.data, cnt.data, ercc.data, zero.data, 
       err.data, batch.id, samp.dat, 
       file = out.fn)
}

