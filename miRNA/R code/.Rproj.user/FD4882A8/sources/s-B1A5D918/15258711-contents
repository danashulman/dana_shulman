#countsnorm <- one quant file
library( "biomaRt" )
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",version = "Ensembl Genes 100", host = "uswest.ensembl.org")
#listAttributes(ensembl)
#listEnsembl()

#some times doesn't work because the server is full, so switch host (^) and try again.
genemap1 <- getBM(attributes = c("ensembl_transcript_id_version","ensembl_transcript_id" ,"ensembl_gene_id", "ensembl_gene_id_version","entrezgene_id","external_gene_name",
                                 "gene_biotype","description","chromosome_name","start_position","end_position"),
                  filters = "ensembl_transcript_id_version",
                  values = rownames(countsnorm),
                  mart = ensembl)

idx_counts <- match(rownames(countsnorm), genemap1$ensembl_transcript_id_version )

countsnorm_DF <- as.data.frame(countsnorm[, 1])

countsnorm_DF$ensembl_gene_id <- genemap1$ensembl_gene_id[idx_counts]
countsnorm_DF$ensembl_gene_id_version <- genemap1$ensembl_gene_id_version[idx_counts]
countsnorm_DF$entrez <- genemap1$entrezgene_id[idx_counts]
countsnorm_DF$external_gene_name <- genemap1$external_gene_name[idx_counts]
countsnorm_DF$rnacentral <- genemap1$rnacentral[idx_counts]
countsnorm_DF$gene_biotype <- genemap1$gene_biotype[idx_counts]
countsnorm_DF$description <- genemap1$description[idx_counts]
countsnorm_DF$chromosome_name <- genemap1$chromosome_name[idx_counts]
countsnorm_DF$start_position <- genemap1$start_position[idx_counts]
countsnorm_DF$end_position <- genemap1$end_position[idx_counts]

write.csv(countsnorm_DF, "list_gene.csv")