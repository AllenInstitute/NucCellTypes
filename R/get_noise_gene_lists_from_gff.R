gff.df <- read.table(file="//allen/programs/celltypes/workgroups/rnaseqanalysis/RNAseq/indexes/GRCm38/20160118_GRCm38.p3.gtf", sep = "\t", stringsAsFactors = FALSE)
gff.df$gene <- sapply(gff.df$V9, function(x) {
  gene1 <- strsplit(as.character(x), ";")[[1]][2]
  sub(" gene_symbol ", "", gene1)
})
gff.df$gene_id <- sapply(gff.df$V9, function(x) {
  gene1 <- strsplit(as.character(x), ";")[[1]][1]
  sub(" gene_id ", "", gene1)
})

# X/Y chr genes
sex.genes.df <- subset(gff.df, V1 %in% c("NC_000086.7", "NC_000087.7"))  # X, Y chr
sex.genes <- unique(sex.genes.df$gene)
write.table(sex.genes, file = "//allen/programs/celltypes/workgroups/humancelltypes/Analyses/CellTypes/Cell_vs_Nuc/lib/sex_genes_20160118_GRCm38.p3.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Mito encoded genes
mito.genes.df <- subset(gff.df, V1 %in% c("NC_005089.1"))  # mito genome
mito.genes <- unique(mito.genes.df$gene)
write.csv(mito.genes, file = "//allen/programs/celltypes/workgroups/humancelltypes/Analyses/CellTypes/Cell_vs_Nuc/lib/mito_genes_20160118_GRCm38.p3.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
