FixDateGenes <- function(x) {
  x <- as.character(x)
  x.num <- gsub("-.*", "", x)
  x.new <- ifelse(grepl("-Mar$", x), paste0("March", x.num), 
                     ifelse(grepl("-Sep$", x), paste0("Sept", x.num), 
                            x))
  x.new2 <- ifelse(x.new == "Sept15", "Sep15", x.new)
}

gene.df <- read.table(file = "lib/mouse_GRCm38_start_stop_SHORT.DAT", sep = "\t")
colnames(gene.df) <- c("gene_id", "start_bp", "end_bp", "V4", "gene")
gene.df$gene <- sapply(gene.df$gene, FixDateGenes)
gene.df$gene_len <- gene.df$end_bp - gene.df$start_bp + 1
load("lib/exonic.gene.sizes.Rdata")
transcript.len <- unlist(exonic.gene.sizes)
gene.df$transcript_len <- transcript.len[match(gene.df$gene, names(transcript.len))]

fix.genes <- which(gene.df$gene_len < gene.df$transcript_len)
gene.df$gene_len[fix.genes] <- gene.df$transcript_len[fix.genes]

gene.df$intron_len <- gene.df$gene_len - gene.df$transcript_len

gene.df.unique <- unique(gene.df[, c("gene", "gene_len", 
                                     "transcript_len", "intron_len")])

write.csv(gene.df.unique,
          file = "lib/mouse_GRCm38_gene_len.csv", row.names = FALSE)
