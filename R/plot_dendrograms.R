load("C:/Users/trygveb/Dropbox/AIBS/Transcriptomics/Manuscripts/WholeCell_vs_Nuc/Figures/Figure 4/dendrograms.rdata", verbose = TRUE)

pdf(file = "C:/Users/trygveb/Dropbox/AIBS/Transcriptomics/Manuscripts/WholeCell_vs_Nuc/Figures/Figure 4/dendrograms.pdf", width = 7, height = 4)
plot(dendI, las = 1, ylim = c(0, 1))
plot(dendE, las = 1, ylim = c(0, 1))
dev.off()