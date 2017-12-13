root.path <- "//allen/programs/celltypes/workgroups"

# Load functions
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(dendextend)

source(file = paste0(root.path, "/rnaseqanalysis/Script_Repository/Trygve/src/iterclust/fCalcCC.R"))
source(file = paste0(root.path, "/rnaseqanalysis/Script_Repository/Trygve/src/iterclust/fCutCC.R"))
source(file = paste0(root.path, "/rnaseqanalysis/Script_Repository/Trygve/src/fMergeSimilarClusters.R"))
source(file = paste0(root.path, "/rnaseqanalysis/Script_Repository/Trygve/src/fMakeClNetwork.R"))


for (samp.set in c("nuc", "cell", "snap25")[1]) {
  for (read.set in c("intron_exon", "exon", "intron", "varIE_clE", "varE_clIE")[1]) {
    project.name <- "mouse_cell_vs_nuc"
    dat.path <- paste0(root.path, "/hct/CT_clustering/", project.name)
    batch.ext <- ifelse(read.set == "intron_exon", samp.set, 
                        paste0(samp.set, "_", read.set))
    batch.date <- paste0("20171211_VISp_", batch.ext)
    print(batch.date)
    run.num <- 1
    analysis.path <- paste0(dat.path, "/cache/", batch.date, "/", run.num)
    path1 <- paste0(dat.path, "/cache/", batch.date, "/", run.num, "/")
    
    # Generate and cluster co-clustering matrix
    CalcCC(pl.path = path1, skip.dir = NULL)
    CutCC(root.path = root.path, cl.cons.path = path1, gap.seq = seq(0.01, 0.99, 0.01))
    
    # Select gap that gives median number of clusters
    cl.num <- read.csv(paste0(path1, "cl_num.csv"))
    cl.median <- median(cl.num$x)
    cc.pw.diff <- read.csv(paste0(path1, "cc.pw.diff.csv"))
    cut.height <- max(cc.pw.diff$gap1[which(cc.pw.diff$num.cl - cl.median >= 0)])
    CutCC(root.path = root.path, cl.cons.path = path1, gap.seq = cut.height)
    
    
    #### Curate clusters ####
    
    #### Load data ####
    load(paste0(dat.path, "/data/", batch.date, "_iter_cl_data.rda"), verbose = TRUE)
    cl.cons <- read.csv(paste0(analysis.path, "/cl.cons.csv.gz"), row.names = 1)
    samp.dat$roi[is.na(samp.dat$roi)] <- "VISp"
    samp.dat$roi <- as.factor(samp.dat$roi)
    rownames(nbt.data) <- toupper(rownames(nbt.data))
    
    for (run1 in 1:2) {
      
      #### Add cluster info ####
      cl.path <- paste0(root.path, "/hct/CT_clustering/", project.name, "/cache/", 
                        batch.date, "/", run.num, "/cl_merge.csv")
      cl.bkp.path <- paste0(root.path, "/hct/CT_clustering/", project.name, "/cache/", 
                            batch.date, "/", run.num, "/cl_merge - orig.csv")
      file.copy(cl.path, cl.bkp.path)
      cl.df <- read.csv(file = cl.path, row.names = 1)
      cl.df$cellid <- sapply(rownames(cl.df), function(x) strsplit(x, "~")[[1]][1])
      samp.dat$exp_component_name[is.na(samp.dat$exp_component_name)] <- cl.df$cellid[is.na(samp.dat$exp_component_name)]
      row.names(samp.dat) <- samp.dat$exp_component_name
      samp.dat$cluster <- as.factor(paste0("cl", cl.df$x[match(samp.dat$exp_component_name, cl.df$cellid)]))
      
      
      
      #### Recompute mean co-clustering matrix ####
      cl.cons.mean <- apply(cl.cons, 1, function(x) tapply(x, samp.dat$cluster, mean))
      cl.cons.mean <- apply(cl.cons.mean, 1, function(x) tapply(x, samp.dat$cluster, mean))
      write.csv(cl.cons.mean, file = paste0(analysis.path, "/cl_merge_robustness.csv"), 
                row.names = TRUE)
      
      within.cc <- diag(cl.cons.mean)
      diag(cl.cons.mean) <- 0
      max.between.cc <- apply(cl.cons.mean, 1, max)
      diag(cl.cons.mean) <- within.cc
      
      # Plot within vs. between co-clustering
      hm.colors <- colorRampPalette(c("white", brewer.pal(9, "YlOrRd")))(100)
      
      pdf(file = paste0(analysis.path, "/cl_merge_robustness.pdf"))
      pheatmap(cl.cons.mean, color = hm.colors, fontsize = 4, border_color = NA)
      plot(max.between.cc, within.cc, cex = 0.3 * sqrt(table(samp.dat$cluster)), 
           xlim = c(0, 1), ylim = c(0, 1), 
           xlab = "Max between cluster co-clustering", 
           ylab = "Average within cluster co-clustering",
           main = "Cluster robustness")
      points(median(max.between.cc), median(within.cc), pch = 19, col = "red", cex = 2)
      abline(0, 1)
      abline(0.2, 1, lty = "dotted")
      dev.off()
      
      
      #### Outlier clusters ####
      # Flag clusters based on marker gene expression
      genes <- c("SNAP25", "GAD1", "GAD2", "SLC17A7", "SLC17A6")
      rownames(nbt.data) <- toupper(rownames(nbt.data))
      expr.cnt <- apply(nbt.data[genes, ], 1, function(x) {
        expr.thresh <- 1
        tapply(x, samp.dat$cluster, function(x) sum(x > expr.thresh) / length(x))
      })
      
      pdf(file = paste0(analysis.path, "/cluster_marker_check.pdf"), 
          width = 8, height = 8)
      pairs(expr.cnt, main = "Expression QC check")
      dev.off()
      
      expr.outlier <- (expr.cnt[, "GAD1"] > 0.5 | expr.cnt[, "GAD2"] > 0.5) & 
        (expr.cnt[, "SLC17A7"] > 0.5 | expr.cnt[, "SLC17A6"] > 0.5)
      
      # Flag clusters based on QC metrics
      neun.pos.cnt <- expr.cnt[, "SNAP25"]
      neun.thresh <- 0.75
      
      qc.metrics <- c("percent_reads_aligned_total",
                      "percent_reads_unique")
      qc.metrics <- intersect(qc.metrics, colnames(samp.dat))
      qc.median <- apply(samp.dat[, qc.metrics], 2, function(x) {
        tapply(x, samp.dat$cluster, function(y) median(as.numeric(y), na.rm = TRUE))
      })
      CheckOutlier <- function(x, iqr.mult = 3) {
        med.diff <- median(x, na.rm=TRUE) - x
        iqr <- abs(quantile(x, 0.75, na.rm=TRUE) - quantile(x, 0.25, na.rm=TRUE))
        outlier.check <- ifelse(med.diff > iqr.mult * iqr, 1, 0)
        return(outlier.check)
      }
      qc.check <- apply(qc.median, 2, CheckOutlier, 3)
      
      if (sd(qc.check) > 0) {
        pdf(file = paste0(analysis.path, "/cluster_qc_check.pdf"),
            width = 8, height = 16, onefile = FALSE)
        pheatmap(qc.check, cluster_rows = FALSE, color = c("grey90", "black"))
        dev.off()
      }
      
      qc.outlier <- neun.pos.cnt > neun.thresh & apply(qc.check, 1, sum) > 0
      
      # Identify outliers
      outlier.cl <- unique(c(names(qc.outlier)[qc.outlier == TRUE], 
                             names(expr.outlier)[expr.outlier == TRUE]))
      
      
      #### Define cluster classes #####
      # Group clusters - exc, inhib, glia, outliers
      cl.class <- list()
      cl.class[["inh"]] <- setdiff(names(neun.pos.cnt)[neun.pos.cnt > neun.thresh & 
                                                         (expr.cnt[, "GAD1"] > 0.5 | 
                                                            expr.cnt[, "GAD2"] > 0.5)], 
                                   c(outlier.cl))
      cl.class[["exc"]] <- setdiff(names(neun.pos.cnt)[neun.pos.cnt > neun.thresh & 
                                                         (expr.cnt[, "SLC17A7"] > 0.2 |
                                                            expr.cnt[, "SLC17A6"] > 0.2)], 
                                   c(outlier.cl, cl.class[["inh"]]))
      cl.class[["glia"]] <- setdiff(names(neun.pos.cnt)[neun.pos.cnt <= neun.thresh], 
                                    c(outlier.cl))
      cl.class[["outlier"]] <- setdiff(levels(samp.dat$cluster), 
                                       unlist(cl.class[c("inh", "exc", "glia")]))
      
      #### Write cluster annotation file ####
      cl.anno <- data.frame(cluster_order = 1:length(unlist(cl.class)), 
                            cluster_id = unlist(cl.class), 
                            cluster_type = gsub("[0-9]*", "", names(unlist(cl.class))))
      cl.anno.path <- paste0(root.path, "/hct/CT_clustering/", project.name, "/cache/", 
                             batch.date, "/", run.num, "/cl.anno.csv")
      write.csv(cl.anno, file = cl.anno.path, row.names = FALSE)
      
      # Co-clustering diagrams
      MakeClNetwork(root.path = root.path, project.name = project.name, 
                    batch.date = batch.date, run.num = run.num)
      
      
      
      
      if (run1 == 1) {
        
        #### Merge clusters ####
        cl.df.new <- cl.df
        gene.sd <- apply(nbt.data, 1, sd)
        
        cl.set.names <- list(c("inh", "exc", "glia"), "outlier")
        cl.sets <- sapply(cl.set.names, function(x) na.omit(unlist(cl.class[x])))
        cl.sets <- cl.sets[sapply(cl.sets, length) > 1]
        
        for (set1 in cl.sets) {
          min.markers <- 1
          branch.height <- 0.5
          cc.diff.thresh <- 0.25
          merge.cnt <- 1
          while (merge.cnt > 0) {
            samp.dat$cluster <- as.factor(paste0("cl", cl.df.new$x[match(rownames(samp.dat), 
                                                                         cl.df.new$cellid)]))
            cl.expr <- t(apply(nbt.data[order(-gene.sd)[1:1000], ], 1, 
                               function(x) tapply(x, samp.dat$cluster, median)))
            cl.subset <- intersect(colnames(cl.expr), set1)
            if (length(cl.subset) <= 1) {
              # Stop if only one cluster left
              merge.cnt <- 0
            } else {
              expr.cor <- cor(cl.expr[, cl.subset])
              cl.tree <- hclust(as.dist(1 - expr.cor))
              cl.groups <- cutree(cl.tree, h = branch.height)
              
              if (all(table(cl.groups) == 1)) {
                # Stop if all clusters above branch height
                merge.cnt <- 0
              } else {
                # Calc DEX gene pairs and plot markers of closest pairs
                min.dex.df <- MergeSimilarClusters(root.path = "//allen", cl.cons, 
                                                   samp.dat, nbt.data, batch.date, cl.groups, 
                                                   mean.expr.thresh = 1)
                
                # Merge clusters with < N robust markers
                merge.cnt <- sum(min.dex.df$min.dex < min.markers | 
                                   min.dex.df$min.cc.diff < cc.diff.thresh)
                print(merge.cnt)  # Track progress
                merged.cl <- NULL
                if (merge.cnt > 0) {
                  cl.tomerge <- droplevels(subset(min.dex.df, min.dex < min.markers |
                                                    min.cc.diff < cc.diff.thresh))
                  for (i in 1:nrow(cl.tomerge)) {
                    fromcl <- cl.tomerge$`2`[i]
                    tocl <- cl.tomerge$`1`[i]
                    if (! (fromcl %in% merged.cl | tocl %in% merged.cl)) {
                      cl.df.new$x[which(cl.df.new$x == fromcl)] <- tocl
                      merged.cl <- c(merged.cl, fromcl, tocl)
                    }
                  }
                }
              }
            } 
          }
        }
        length(table(cl.df.new$x))
        
        # Save cluster tree
        expr.cor <- cor(cl.expr)
        cl.tree <- as.dendrogram(hclust(as.dist(1 - expr.cor)))
        class.color <- list("orange", "blue", "green", "red")
        names(class.color) <- names(cl.class)
        color1 <- unlist(mapply(rep, class.color, sapply(cl.class, length)))
        labels_colors(cl.tree) <- color1[match(colnames(expr.cor)[order.dendrogram(cl.tree)], 
                                               unlist(cl.class))]
        pdf(file = paste0(analysis.path, "/cluster_tree.pdf"), 
            width = 12, height = 4)
        plot(cl.tree, ylab = "Expression distance (1 - r)")
        legend("topright", fill = unlist(class.color), legend = names(class.color), bty = "n")
        dev.off()
        
        
        write.csv(cl.df.new, file = cl.path)
        
      }
    }
  }
}
