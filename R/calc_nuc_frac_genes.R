#### Calc nuclear fraction of all genes ####
load("//allen/programs/celltypes/workgroups/hct/CT_clustering/mouse_cell_vs_nuc/scriptsAndOutputJAM/outputFiles_new/info.RData", verbose = TRUE)
type.nuc.frac <- read.csv("C:/Users/trygveb/Dropbox/AIBS/Transcriptomics/Manuscripts/WholeCell_vs_Nuc/Tables/nuclearFraction_snoRNA.csv")

nuc.ratio <- matrix(NA, nrow = nrow(meansIN), ncol = 2, 
                    dimnames = list(row.names(meansIN), 
                                    c("nucfrac", "nucfrac_sd")))
for (i in 1:nrow(meansIN)) {
  ratio.mean <- meansIN[i, ] / meansIC[i, ] * type.nuc.frac$nucfrac_mean
  excl.types <- which(meansIN[i, ] == 0 | meansIC[i, ] == 0)
  ratio.mean[excl.types] <- NA
  ratio.var <- ratio.mean^2 * ((sdsIN[i, ] / meansIN[i, ])^2 +
                                 (sdsIC[i, ] / meansIC[i, ])^2 +
                                 (type.nuc.frac$nucfrac_mean_sd / 
                                    type.nuc.frac$nucfrac_mean)^2)
  
  ratio.mean[ratio.mean > 1 & (ratio.mean - sqrt(ratio.var)) < 1] <- 1
  # ratio.wts <- apply(cbind(meansIN[i, ], meansIC[i, ]), 1, max)
  ratio.wts <- meansIC[i, ]
  ratio.wts[excl.types] <- 0
  ratio.wts <- ratio.wts / sum(ratio.wts)
  wt.ratio.mean <- weighted.mean(ratio.mean, ratio.wts, na.rm = TRUE)
  wt.ratio.var <- sum(ratio.wts * (ratio.mean - wt.ratio.mean)^2, na.rm = TRUE)
  
  # Check if 95% CI spans range
  wt.ratio.q05 <- wt.ratio.mean - 2 * sqrt(wt.ratio.var)
  wt.ratio.q95 <- wt.ratio.mean + 2 * sqrt(wt.ratio.var)
  excl.ratio <- wt.ratio.mean > 1 | (wt.ratio.q05 < 0 & wt.ratio.q95 > 1)
  if (is.na(excl.ratio) | excl.ratio) {
    nuc.ratio[i, ] <- c(NA, NA)
  } else {
    nuc.ratio[i, ] <- c(wt.ratio.mean, sqrt(wt.ratio.var))
  }
}

