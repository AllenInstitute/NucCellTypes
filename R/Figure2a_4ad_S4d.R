###################################################################################################################
# This script needs to be run from a windows computer with write-access in R version 3.0 or higher.

print("Set working folders, read in data, and load libraries.")

mainFolder     = ""
outputFolder   = paste(mainFolder,"output/",sep="")
scriptsFolder  = paste(mainFolder,"R/",sep="")
inputFolder    = paste(mainFolder,"data/",sep="")

# Load these libraries
options(stringsAsFactors=FALSE)

library(beeswarm)
library(WGCNA)
library(edgeR)   
library(feather)
library(dendextend)
library(monocle)
library(ggplot2)
library(dplyr)
library(matrixStats)
library(subSeq)

# Read in the data and extra scripts
source(paste(scriptsFolder,"extraFunctions.r",sep=""))
load(paste0(inputFolder, "start_data.rda"))
 # This includes "annoCel", "annoNuc", "countsR", "erccTable", "exprAll", "introns", "sampAll"
 # the "All" samples include all cell and nuclei, and 30 "ContronTotalRNA" samples
cpmE = cpm(countsR)
cpmI = cpm(introns)


#########################################################################################
## FIGURE 2A ----------------------------------------------------------------------------
## FIGURE 2A ----------------------------------------------------------------------------

print("Figure 2A: Plot overall statistics for nuclei, cells, and controls.")

kpAlignStats   = c("percent_reads_aligned_to_exons","percent_reads_aligned_to_rrna","percent_reads_aligned_to_trna","percent_reads_aligned_to_other_ncrna","percent_reads_aligned_to_mt_exons","percent_reads_aligned_to_mt_rrna","percent_reads_aligned_to_mt_trna","percent_reads_aligned_to_mt_other_ncrna","percent_reads_aligned_intron","percent_reads_intergenic","percent_reads_aligned_to_ecoli","percent_reads_aligned_to_synthetic_constructs")
 
otherVars      = c("cre_line","batch_vendor_name","Fill.Date",kpAlignStats,"ERCCLimit")

sampleType  = sampAll$Type
alignStatsA = sampAll[,c(kpAlignStats,"percent_unmapped_reads")]
kpAlignStats2 = c("percent_unmapped_reads","percent_reads_aligned_to_exons","percent_reads_aligned_intron",
    "percent_reads_aligned_to_other_ncrna","percent_reads_intergenic")
alignStats = alignStatsA[,kpAlignStats2]
alignStats$percent_reads_aligned_to_everything_else = 100-rowSums(alignStats,na.rm=TRUE)

nonGenome   = c("percent_reads_aligned_to_rrna","percent_reads_aligned_to_trna","percent_reads_aligned_to_mt_exons",
     "percent_reads_aligned_to_mt_rrna","percent_reads_aligned_to_mt_trna","percent_reads_aligned_to_mt_other_ncrna",
	 "percent_reads_aligned_to_ecoli","percent_reads_aligned_to_synthetic_constructs")

genomePerc  = sampAll$percent_reads_aligned_total-rowSums(sampAll[,nonGenome])
genomeReads = round(sampAll$total_reads*(1/100)*genomePerc)

EII = data.frame(ExonReads=round(colSums(countsR)),IntronReads = colSums(introns), 
  IntergenicReads = genomeReads-colSums(countsR)-colSums(introns))
readLocs   = 100*EII/rowSums(EII)
readLocs$PercentAlignedToGenome = genomePerc
readLocs = readLocs[,c(4,1:3)]
rownames(readLocs) <- colnames(exprAll)

Ns = c(0,1,2,5,10,25,50,100,250,500,1000)
geneCounts = list()
for (n in Ns){
 N = as.character(n)
 geneCounts[[N]] = data.frame(AnyCounts = colSums((cpmE>n)|(cpmI>n)), 
   ExonCounts = colSums(cpmE>n), IntronCounts = colSums(cpmI>n))
 rownames(geneCounts[[N]]) <- colnames(exprAll)
}
# For the manuscript, CPM > 0 is shown.

stats = geneCounts
stats[["ALIGNMENT_STATS"]] = readLocs
cols = c("#387EB8","#E21F26","grey")
yMax = c(rep(13000,6),rep(6000,2),rep(2000,3),100)
names(yMax) = names(stats)

pdf(paste0(outputFolder, "plotAllStats_beeswarm.pdf"),height=10,width=4.5)
for (N in names(stats)){
 for (cl in colnames(stats[[N]])){
  out = data.frame(values = stats[[N]][,cl],SampleType = as.factor(sampleType))
  beeswarm(values ~ SampleType, data = out, 
    col = cols, pch = 16, corral = "random", cex=0.9, method="center",
    main = paste(cl,N,sep="-"),ylim = c(0,yMax[N]),las=1,cex.axis=1.35,ylab="")
  bxplot(values ~ SampleType, data = out, add = TRUE)
  abline(h=0)
 }
}
dev.off()
# Pages 1-3 and 34-37 are in the manuscript

print("===== Determine statistical significance.")

pvalAnova <- cbind(readLocs,geneCounts[[1]])
pvalsAll  <- apply(pvalAnova,2,getAnovaPvalforApply,sampleType)
norc      <- is.element(sampleType,c("Cells","Nuclei"))
pvalsNorc <- apply(pvalAnova[norc,],2,getAnovaPvalforApply,sampleType[norc])
cbind(pvalsAll,pvalsNorc)
#                             pvalsAll     pvalsNorc
# PercentAlignedToGenome  5.708036e-41  8.199823e-14
# ExonReads               3.758943e-95  4.077768e-87
# IntronReads            9.712306e-127 3.175153e-105
# IntergenicReads        3.541679e-178 2.043684e-112
# AnyCounts               0.000000e+00  0.000000e+00
# ExonCounts              0.000000e+00  0.000000e+00
# IntronCounts           2.794476e-148 6.719808e-109


#########################################################################################
## FIGURE 4 -----------------------------------------------------------------------------
## FIGURE 4 -----------------------------------------------------------------------------

print("Figure 4A/D: Build clustering trees for cells, nuclei, introns, and exons.")

kpCell    = sampleType=="Cells"
kpNuc     = sampleType=="Nuclei"
kpBoth    = kpCell|kpNuc
cl        = paste(substr(sampleType,1,1),sampAll$cluster,sep="_")
names(cl) = rownames(sampAll)

# For ordering the dendrograms as in the paper
clusters= sort(unique(cl))[c(10:11,1:9,21:22,20,12:19)]  # Clusters, ordered by layer
l.rank  = setNames(1:22, clusters)
cluI    = c("C_Pvalb.Wt1","N_Pvalb.Wt1","C_Sst.Cbln4","N_Sst.Cbln4",
   "C_L6a.Plcxd3","N_L6a.Plcxd3","N_L5b.Cdh13","N_L5b.Samd3",
   "C_L5b.Cdh13","C_L5b.Samd3","C_L5.Chrna6","N_L5.Chrna6",
   "C_L4.Arf5","N_L4.Arf5","C_L4.Hsd11b1","N_L4.Hsd11b1",
   "C_L5a.Batf3","N_L5a.Batf3","C_L6a.Mgp","N_L6a.Mgp",
   "C_L6.Car12","N_L6.Car12")
cluE    = c("C_Pvalb.Wt1","C_Sst.Cbln4","C_L4.Arf5","C_L4.Hsd11b1",
   "C_L5.Chrna6","C_L5b.Cdh13","C_L5b.Samd3","C_L6a.Plcxd3",
   "C_L6a.Mgp","C_L5a.Batf3","C_L6.Car12","N_Pvalb.Wt1",
   "N_Sst.Cbln4","N_L6a.Plcxd3","N_L5.Chrna6","N_L4.Arf5",
   "N_L4.Hsd11b1","N_L5a.Batf3","N_L6.Car12","N_L6a.Mgp",
   "N_L5b.Cdh13","N_L5b.Samd3")

# Lists for the four dendrograms 
clL       <- list(cl[kpNuc],cl[kpCell],cl[kpBoth],cl[kpBoth])
exprDataL <- list(cpm(countsR[,kpNuc]+introns[,kpNuc]),cpm(countsR[,kpCell]+introns[,kpCell]),cpmI[,kpBoth],cpmE[,kpBoth])
lRankL    <- list(l.rank[12:22]-11,l.rank[1:11],setNames(1:22, cluI),setNames(1:22, cluE))
mains     <- c("Nuclei","Cells","Introns","Exons")

# Plot the dendrograms
pdf(paste0(outputFolder, "clusterDendrograms.pdf"))
for (i in 1:length(mains))
   buildAndPlotTree(exprDataL[[i]],clL[[i]],lRankL[[i]],topNgenes=1200,main=mains[i])
dev.off()
# Note that some additional formatting was done in Illustrator (e.g., colored dots, higher-weight lines).


#########################################################################################
## FIGURE S4 -----------------------------------------------------------------------------
## FIGURE S4 -----------------------------------------------------------------------------

print("Figure S4D: Compare intron and exon expression in cells and nuclei.")

pdf(paste0(outputFolder, "intronExonCorrelation_cellsVsNuclei.pdf"))
ne = log2(rowMeans(cpmE[,kpNuc])+1)
ce = log2(rowMeans(cpmE[,kpCell])+1)
verboseScatterplot(ne,ce,xlab="Mean log2_CPM(counts in nuclei)",ylab="Mean log2_CPM(counts in cells)",main="Exons",pch=19,cex=0.5)
ne = log2(rowMeans(cpmI[,kpNuc])+1)
ce = log2(rowMeans(cpmI[,kpCell])+1)
verboseScatterplot(ne,ce,xlab="Mean log2_CPM(counts in nuclei)",ylab="Mean log2_CPM(counts in cells)",main="Introns",pch=19,cex=0.5)
dev.off()


#########################################################################################
## TBD ----------------------------------------------------------------------------
## TBD ----------------------------------------------------------------------------

print("Subsample introns and exons and count the number of genes detected (with counts/CPM>0) at various cutoffs in all cases.")

# Define some variables
sampleType  <- sampAll$Type
counts      <- countsR+introns # TOTAL COUNTS
subCountsN  <- list(All=counts[,sampleType=="Nuclei"])
subCountsC  <- list(All=counts[,sampleType=="Cells"])
medC        <- median(colSums(counts[,sampleType=="Cells"]))
medN        <- median(colSums(counts[,sampleType=="Nuclei"]))
nucScale    <- medC / medN
propNamesC  <- round(medC/10000)/100 
propNamesN  <- round(medN/10000)/100

# Subsample to multiple count levels
proportions <- (11:1)/10
for (i in 2:length(proportions)){
  subCountsN[[i]] <- generateSubsampledMatrix(subCountsN[[1]], proportions[i]*nucScale, seed=i)
  subCountsC[[i]] <- generateSubsampledMatrix(subCountsC[[1]], proportions[i], seed=i)
  propNamesC      <- c(propNamesC,round(medC*proportions[i]/10000)/100)
  propNamesN      <- c(propNamesN,round(medC*proportions[i]/10000)/100)
}

# Calculate genes detected at various cutoffs for cells and nuclei
Ns     <- 0
outNuc <- outCell <- NULL
for (i in 1:length(proportions)){
  outNuc  <- rbind(outNuc,cbind(colSums(subCountsN[[i]]>Ns)/1000,propNamesN[i],"Nuclei"))
  outCell <- rbind(outCell,cbind(colSums(subCountsC[[i]]>Ns)/1000,propNamesN[i],"Cells"))
}
colnames(outNuc) <- colnames(outCell) <- c("GenesDetected","TotalReads","Source")
out <- rbind(outNuc,outCell)
out <- as.data.frame(out)
out$GenesDetected <- as.numeric(out$GenesDetected)
out$TotalReads[out$TotalReads=="2.18"] = "All"
out$TotalReads <- as.factor(as.character(out$TotalReads))
out$Source <- as.factor(out$Source)
rownames(out) <- NULL

# Plot and output the results
pdf(paste0(outputFolder,"genesDetected_subsample.pdf"),height=10,width=8)
outPlot <- ggplot(out, aes(x = TotalReads, y = GenesDetected, fill = Source)) + 
           geom_violin(draw_quantiles=c(0.25,0.5,0.75)) +
		   labs(x="Median number of reads per cell (millions)", y="Genes Detected with CPM>0 (thousands)") +
		   scale_fill_manual(values = c("#387EB8","#E21F26")) +
           scale_y_continuous(breaks=(0:7)*2,limits=c(0,14)) 
outPlot
dev.off()

