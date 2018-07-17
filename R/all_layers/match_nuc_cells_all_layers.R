library(feather)
library(WGCNA)

input.path.cell <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/STARforLIMS/Mouse/Final_Output/Mouse_Final_21913"
input.path.nuc <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/STARforLIMS/Mouse/Final_Output/Nuclei_Noisemodel/R_objects/"
out.path <- "//allen/programs/celltypes/workgroups/hct/CT_clustering/mouse_cell_vs_nuc/cache/20180205_VISp_nuc"

#### Load V1 L5 whole cell data ####
for (file1 in c("samp.dat", "intron", "exon", "zero.wt")) { 
  load(paste0(input.path.cell, "/", file1, "_mouse_star_21913.Rdata"), 
       verbose = TRUE)
}

keep.cell <- which(grepl("VISp", samp.dat$roi))

samp.cell <- droplevels(samp.dat[keep.cell, ])
cd.cell <- introns[, keep.cell] + countsR[, keep.cell]
cpm.cell <- round(sweep(cd.cell, 2, colSums(cd.cell), "/") * 1e6, 0)
expr.cell <- log2(cpm.cell + 1)
wt.cell <- zero.wt$matw[, row.names(samp.cell)]  # Dropout weight matrix calculated using scde as described in Methods


#### Load V1 L5 nuclei data ####
for (file1 in c("samp.dat", "intron", "exon", "zero.wt")) { 
  load(paste0(input.path.nuc, "/", file1, "_MouseNuclei6960_star.Rdata"),
       verbose = TRUE)
}

# Cluster anno
cl.df <- read.csv(file = "//allen/programs/celltypes/workgroups/hct/CT_clustering/mouse_cell_vs_nuc/cache/20180205_VISp_nuc/1/9999/merged/clid_merged.csv", stringsAsFactors = FALSE)
nuc.id <- sapply(cl.df$sampid, function(x) strsplit(x, "~")[[1]][1])

# QC criteria
keep.nuc <- which(row.names(samp.dat) %in% nuc.id)
samp.nuc <- droplevels(samp.dat[keep.nuc, ])
cd.nuc <- intron[, keep.nuc] + exon[, keep.nuc]
cpm.nuc <- round(sweep(cd.nuc, 2, colSums(cd.nuc), "/") * 1e6, 0)
expr.nuc <- log2(cpm.nuc + 1)
wt.nuc <- zero.wt$matw[, row.names(samp.nuc)]

# Keep genes
var.genes.nuc <- apply(expr.nuc, 1, function(x) sum(x > 1) >= 4)
var.genes.cell <- apply(expr.nuc, 1, function(x) sum(x > 1) >= 4)
keep.genes <- intersect(row.names(expr.nuc)[var.genes.nuc], 
                        row.names(expr.cell)[var.genes.cell])
expr.cell <- expr.cell[keep.genes, ]
wt.cell <- wt.cell[keep.genes, ]
expr.nuc <- expr.nuc[keep.genes, ]
wt.nuc <- wt.nuc[keep.genes, ]


#### Find nuc/cell correlations ####
expr.nuc.cell <- cbind(expr.nuc, expr.cell)
cor.nuc.cell <- WGCNA::cor(expr.nuc.cell)

nnuc <- ncol(expr.nuc)
ncell <- ncol(expr.cell)
cor.all <- cor.nuc.cell[1:nnuc, (nnuc + 1):(nnuc + ncell)]
cor.nuc <- cor.nuc.cell[1:nnuc, 1:nnuc]
diag(cor.nuc) <- 0


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


# Select best matching nuclei
cor.nuc.ordered <- cor.nuc[order.cell, order.cell]
max.nuc.all <- apply(cor.nuc.ordered, 1, which.max)
max.cor.nuc <- apply(cor.nuc.ordered, 1, max)


#### Find matching cell/cell correlations ####
cor.cell <- cor.nuc.cell[names(max.cell.all), names(max.cell.all)]
diag(cor.cell) <- 0


# Sample ids
samp.ids <- list()
samp.ids[["nuc"]] <- samp.nuc$exp_component_name[order.cell]
samp.ids[["cell"]] <- names(max.cell.all)
samp.ids[["nuc_matched"]] <- samp.nuc$exp_component_name[max.nuc.all][order.cell]


for (id.set in c("cell", "nuc", "nuc_matched")) {
  write.table(samp.ids[[id.set]], file = paste0(out.path, "/", id.set, "_ids.txt"),
              row.names = FALSE, col.names = FALSE)
}

if (! file.exists(paste0(out.path, "/map_nuc_to_cells.rda"))) {
  save(cor.all.ordered, max.cell.all, max.cor.all, 
       cor.nuc.ordered, max.nuc.all, max.cor.nuc,
       samp.ids, cor.cell,
       file = paste0(out.path, "/map_nuc_to_cells.rda"))
}
