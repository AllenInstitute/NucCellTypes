library(beeswarm)
library(feather)
library(dplyr)
library(scrattch)

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/facs/R_Object/20180710_RSC-11-161_mouse_star_samp.dat.Rdata")

rownames(samp.dat)=samp.dat$exp_component_name
samp.dat <- samp.dat[match(sampAll$exp_component_name, samp.dat$exp_component_name), ]
samp.dat$Type <- factor(samp.dat$Type, levels = c("Cells", "Nuclei", "ControlTotalRNA"))


genomic.cols <- c("percent_reads_aligned_total", 
                  "percent_reads_aligned_unique", 
                  "percent_reads_aligned_to_exons", 
                  "percent_reads_aligned_to_introns", 
                  "percent_reads_aligned_to_intergenic") 
# "percent_reads_aligned_to_rrna", 
# "percent_reads_aligned_to_trna", 
# "percent_reads_aligned_to_ncrna"


reads <- matrix(NA, nrow(samp.dat), length(genomic.cols),
                dimnames = list(row.names(samp.dat), genomic.cols))
for (col1 in genomic.cols) {
  reads[, col1] <- samp.dat[, col1] * samp.dat$total_reads
}

reads.pct <- apply(reads, 2, function(x) x / reads[, 2] * 100)
reads.pct[, 1] <- samp.dat$percent_reads_aligned_total
reads.pct[, 2] <- samp.dat$percent_reads_aligned_unique


boxplot(reads.pct[, 1] ~ samp.dat$Type)


cols = c("#387EB8","#E21F26","grey")
# yMax = setNames(c(100, 100, 100, 100), genomic.cols)

pdf("output/plotAllStats_beeswarm.pdf", height=10, width=4.5)
for (col1 in genomic.cols) {
  beeswarm(reads.pct[, col1] ~ samp.dat$Type, 
           corral = "random", method="swarm",
           col = cols, pch = 16, cex=0.9, 
           main = col1, ylim = c(0, 100), 
           las=1, cex.axis=1.35, xlab = "", ylab="")
  bxplot(reads.pct[, col1] ~ samp.dat$Type, add = TRUE)
}
dev.off()





anno$percent_unique_reads_label=samp.dat[anno$sample_id,"percent_reads_aligned_unique"]
unique=anno$percent_unique_reads_label*anno$total_reads_label

anno$percent_reads_aligned_intron_label=samp.dat[anno$sample_id,"percent_reads_aligned_to_introns"]
intron=anno$percent_reads_aligned_intron_label*anno$total_reads_label

anno$percent_reads_aligned_exon_label=samp.dat[anno$sample_id,"percent_reads_aligned_to_exons"]
exon=anno$percent_reads_aligned_to_exons_label*anno$total_reads_label  

intergenic=unique -(exon+intron))  



anno=anno[,-c(grep("intergenic",names(anno)))]

anno$percent_reads_intergenic= 100*(unique -(exon+intron))/unique



anno<-anno%>%annotate_num("percent_reads_intergenic")

#write_feather(anno,"//allen/programs/celltypes/workgroups/hct/CT_clustering/mouse_cell_vs_nuc/shinydb/20170818_VISp_L5_nuc_1/anno_intergenic_modified.feather")

write_feather(anno,"//allen/programs/celltypes/workgroups/hct/CT_clustering/mouse_cell_vs_nuc/shinydb/20170818_VISp_L5_cell_1/anno_intergenic_modified.feather")
