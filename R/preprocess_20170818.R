root.path <- "//allen/programs/celltypes/workgroups"
input.path.cell <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/STARforLIMS/Mouse/Final_Output/Mouse_Final_21913"
input.path.cell.nuc <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/STARforLIMS/Mouse/Final_Output/VISp_L5"
analysis.path <- "//allen/programs/celltypes/workgroups/humancelltypes/Analyses/CellTypes/Cell_vs_Nuc/analysis/20170818/match_nuc_cell"

# Load functions
source('//allen/programs/celltypes/workgroups/rnaseqanalysis/Script_Repository/Trygve/src/fMatchRowCol.R')
source("//allen/programs/celltypes/workgroups/rnaseqanalysis/Script_Repository/Trygve/src/celltypes/fMaptoV1Tree.R")


# Noise genes
sex.genes <- scan(file = "//allen/programs/celltypes/workgroups/humancelltypes/Analyses/CellTypes/Cell_vs_Nuc/lib/sex_genes_20160118_GRCm38.p3.txt", "character")
mito.genes <- scan(file = "//allen/programs/celltypes/workgroups/humancelltypes/Analyses/CellTypes/Cell_vs_Nuc/lib/mouse_mito_genes.txt", "character")

for (samp.set in c("nuc", "cell", "snap25")[2]) {
  for (read.set in c("intron_exon", "exon", "intron", "varIE_clE", "varE_clIE")[4:5]) {
    
    # Load data
    if (samp.set == "nuc") {
      for (file1 in c("samp.dat", "introns", "exon", "zero.wt")) { 
        load(paste0(input.path.cell.nuc, "/Nuc/", file1, "_Mouse_Star_VISp_L5_Nuc.Rdata"))
      }
    } else {
      for (file1 in c("samp.dat", "intron", "exon")) { 
        load(paste0(input.path.cell, "/", file1, "_mouse_star_21913.Rdata"))
      }
      load(paste0(input.path.cell, "/tmp/zero.wt_mouse_star_21913.Rdata"))
    }
    
    # Remove sex and mito genes
    keep.genes <- which(! rownames(introns) %in% c(sex.genes, mito.genes))
    introns <- introns[keep.genes, ]
    countsR <- countsR[keep.genes, ]
    
    
    # Keep QC passed cells
    samp.ids <- scan(file = paste0(analysis.path, "/", samp.set, "_ids.txt"), 
                     "character")
    qc.samp <- which(samp.dat$exp_component_name %in% samp.ids)
    
    # Keep samples/genes with information in all matrices
    samp.dat <- droplevels(samp.dat[qc.samp, ])
    
    # Create expression matrix
    if (read.set == "intron_exon") {
      cd <- introns[, qc.samp] + countsR[, qc.samp]
      cd2 <- cd
    } else if (read.set == "exon") {
      cd <- countsR[, qc.samp]  
      cd2 <- cd
    } else if (read.set == "intron") {
      cd <- introns[, qc.samp] 
      cd2 <- cd
    } else if (read.set == "varIE_clE") {
      cd <- introns[, qc.samp] + countsR[, qc.samp]
      cd2 <- countsR[, qc.samp]  
    } else if (read.set == "varE_clIE") {
      cd <- countsR[, qc.samp]  
      cd2 <- introns[, qc.samp] + countsR[, qc.samp]
    }
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
    samp.fn <- ifelse(read.set == "intron_exon", samp.set, paste0(samp.set, "_", read.set))
    out.fn <- paste0(dat.path, "/20170818_VISp_L5_", samp.fn, "_iter_cl_data.rda")
    if (! file.exists(out.fn)) {
      save(nbt.data, cnt.data, ercc.data, zero.data, 
           err.data, batch.id, samp.dat, 
           file = out.fn)
    }
    
  }
}