root.path <- "//allen/programs/celltypes/workgroups"


library(feather)

source(file = paste0(root.path, "/rnaseqanalysis/Script_Repository/Trygve/src/iterclust/fMakeShinyFeatherDb.R"))

calc_beta <- function(y) {
  d1 <- as.matrix(dist(y))
  eps1 <- 1e-10
  beta1 <- sum(d1^2) / (sum(d1) + eps1)
  return(beta1)
}

calc_tau <- function(vec) {
  xi_vec = vec/max(vec)
  tau = sum(1 - xi_vec) / (length(vec) - 1)
  tau[is.na(tau)] <- 0
  return(tau)
}



for (samp.set in c("nuc", "cell", "snap25")[1]) {
  for (read.set in c("intron_exon", "exon", "intron", "varIE_clE", "varE_clIE")[1]) {
    project.name <- "mouse_cell_vs_nuc"
    dat.path <- paste0(root.path, "/hct/CT_clustering/", project.name)
    batch.ext <- ifelse(read.set == "intron_exon", samp.set, 
                        paste0(samp.set, "_", read.set))
    batch.date <- paste0("20171211_VISp_", batch.ext)
    print(batch.date)
    run.num <- 1
    
    
    # Make Shiny db
    dat.fn <- paste0(dat.path, "/data/", batch.date, "_iter_cl_data.rda")
    MakeShinyFeatherDb(root.path = root.path, dat.path = dat.path, dat.file = dat.fn, 
                       expr.var = "cnt.data", samp.var = "samp.dat", 
                       batch.date = batch.date, run.num = run.num, 
                       db.dir = paste0(dat.path, "/shinydb"), 
                       load.cluster.anno = TRUE, 
                       include.num.anno = TRUE, 
                       anno.cols = c("batch_vendor_name", "gender", 
                                     "facs_population_plan", "cre_line", 
                                     "transgenic_recombinase", "transgenic_reporter", 
                                     "reporter", "roi", "percent_cd_longer_than_400bp", 
                                     "rna_amplification_pass_fail", "library_prep_pass_fail", 
                                     "max.leaf", "max.leaf.class", 
                                     "sample_name", "exp_component_name"),
                       anno.names = NULL)
    
    
    #### Make prop feather ####
    shiny.dir <- paste0(dat.path, "/shinydb/", batch.date, "_", run.num, "/")
    
    Expr.dat <- read_feather(paste0(shiny.dir,"data.feather"))
    Samp.dat <- read_feather(paste0(shiny.dir,"anno.feather"))
    
    datIn2    = as.matrix(Expr.dat[,colnames(Expr.dat)!="sample_id"])
    rownames(datIn2) <- Expr.dat$sample_id
    datIn2    = t(datIn2)
    
    cll        = as.character(Samp.dat$cluster_label)
    names(cll) = Samp.dat$sample_id
    cll        = factor(cll,levels=unique(cll))
    
    exprThresh = 1  # COULD ALSO BE 0
    propmat = do.call("cbind", tapply(names(cll), cll, 
                                      function(x) rowMeans(datIn2[,x]>exprThresh))) 
    write_feather(data.frame(propmat),paste0(shiny.dir,"osnat.prop1.feather"))
    propmat = t(propmat)
    write_feather(data.frame(propmat),paste0(shiny.dir,"prop.feather"))
    
    
    #### Make FPKM feather ####
    input.path.cell <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/STARforLIMS/Mouse/Final_Output/Mouse_Final_21913"
    input.path.cell.nuc <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/STARforLIMS/Mouse/Final_Output/Nuclei_Noisemodel/R_objects"
    
    if (samp.set == "nuc") {
      load(paste0(input.path.cell.nuc, "/fpkm_MouseNuclei6960_star.Rdata"))
    } else {
      load(paste0(input.path.cell, "/fpkm_mouse_star_21913.Rdata"))
    }
    
    fpkmmat = do.call("cbind", tapply(names(cll), cll, function(x) rowMeans(datExprR[,x])))
    write_feather(data.frame(fpkmmat),paste0(shiny.dir,"osnat.log.median.feather"))
    fpkmmat = t(fpkmmat)
    write_feather(data.frame(fpkmmat),paste0(shiny.dir,"fpkm.feather"))
    
    
    #### Calc marker scores ####
    anno <- as.data.frame(Samp.dat)
    keep.cl <- as.character(unique(anno$cluster_label[anno$cluster_type_label %in% 
                                                        c("exc", "inh", "glia")]))
    
    prop1 <- t(as.data.frame(propmat))
    colnames(prop1) <- unique(as.character(anno$cluster_label))
    
    expr1 <- t(as.data.frame(fpkmmat))
    colnames(expr1) <- unique(as.character(anno$cluster_label))
    
    keep.genes <- intersect(rownames(prop1), rownames(expr1))
    prop1 <- prop1[keep.genes, keep.cl]
    expr1 <- expr1[keep.genes, keep.cl]
    
    fpkm.max <- apply(expr1, 1, max)
    fpkm.max.cluster <- colnames(expr1)[apply(expr1, 1, which.max)]
    
    tau.prop1 <- apply(prop1, 1, calc_tau)
    tau.expr <- apply(log2(expr1 + 1), 1, calc_tau)
    beta.prop1 <- apply(prop1, 1, calc_beta)
    beta.expr <- apply(log2(expr1 + 1), 1, calc_beta)
    total.clusters <- apply(prop1, 1, function(x) sum(x > 0.5))
    inh.clusters <- apply(prop1[, as.character(unique(anno$cluster_label[anno$cluster_type_label == "inh"]))], 1, function(x) sum(x > 0.5))
    exc.clusters <- apply(prop1[, as.character(unique(anno$cluster_label[anno$cluster_type_label == "exc"]))], 1, function(x) sum(x > 0.5))
    glia.clusters <- apply(prop1[, as.character(unique(anno$cluster_label[anno$cluster_type_label == "glia"]))], 1, function(x) sum(x > 0.5))
    
    score.df <- data.frame(gene = names(total.clusters), 
                           total.clusters, inh.clusters, exc.clusters, glia.clusters, 
                           fpkm.max.cluster, fpkm.max, 
                           beta.prop1, beta.expr, tau.prop1, tau.expr)
    score.df$gene <- sapply(score.df$gene, function(x) sub("^X([0-9])", "\\1", x))
    write.csv(score.df, file = paste0(shiny.dir, "marker.scores.csv"), row.names = FALSE)
    
  }
}