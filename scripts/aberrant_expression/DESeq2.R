library(tidyverse)
library(stringr)
library(OUTRIDER)
library(DESeq2)

# Load data -------------------------------------------------------------------------------------------

load("/home/dzhang/projects/ATG7_rob_t_analysis/results/get_gene_count_RSE/gene_counts_rse.rda")

ref <- dasper:::.ref_load(ref = "/data/references/ensembl/gtf_gff3/v97/Homo_sapiens.GRCh38.97.gtf")

gene_count_test <- read_delim("results/get_gene_count_RSE/control_1.gene_reads.gct",
                              delim = "\t", skip = 2)

# Main ------------------------------------------------------------------------------------------------

gene_count_ods <- OutriderDataSet(se = gene_counts_rse)

##### Filter for expressed genes #####

gene_count_ods <- filterExpression(gene_count_ods, 
                                   gtfFile = ref, 
                                   filterGenes = FALSE, 
                                   savefpkm = TRUE)

gene_count_ods <- gene_count_ods[mcols(gene_count_ods)$passedFilter,]

##### Run DESeq2 - ATG7 vs rest patients #####

register(MulticoreParam(20))

which_ATG7 <- which(colData(gene_count_ods)[["gene_name"]] == "ATG7")
gene_count_ods$condition <- "control"
gene_count_ods$condition[which_ATG7] <- "case"

gene_count_ods$batch <- as.factor(gene_count_ods$batch)

dds <- DESeqDataSetFromMatrix(countData = assays(gene_count_ods)[["counts"]],
                              colData = colData(gene_count_ods),
                              design = ~ batch + sex + condition)

# run deseq with default settings 
dds <- DESeq(dds, parallel = TRUE)
res <- results(dds, parallel = TRUE)
res <- res[order(res$padj),]

res <- res %>% 
  as.data.frame() %>% 
  rownames_to_column("gene_id") %>%
  left_join(gene_count_test, by = c("gene_id" = "Name")) %>% 
  as_tibble()

ATG7_deseq2_res <- res

save(ATG7_deseq2_res, file = "results/DESeq2/ATG7_deseq2_res.rda")

##### Run DESeq2 - control vs rest patients #####

register(MulticoreParam(20))

which_COL6A <- which(str_detect(colData(gene_count_ods)[["gene_name"]], "COL6A*"))
gene_count_ods$condition <- "control"
gene_count_ods$condition[which_COL6A] <- "case"
  
gene_count_ods$batch <- as.factor(gene_count_ods$batch)

dds <- DESeqDataSetFromMatrix(countData = assays(gene_count_ods)[["counts"]],
                              colData = colData(gene_count_ods),
                              design = ~ batch + sex + condition)

# run deseq with default settings 
dds <- DESeq(dds, parallel = TRUE)
res <- results(dds, parallel = TRUE)
res <- res[order(res$padj),]

res <- res %>% 
  as.data.frame() %>% 
  rownames_to_column("gene_id") %>%
  left_join(gene_count_test, by = c("gene_id" = "Name")) %>% 
  as_tibble()

COL6A_deseq2_res <- res

save(COL6A_deseq2_res, file = "results/DESeq2/COL6A_deseq2_res.rda")

