library(tidyverse)
library(stringr)
library(OUTRIDER)
library(DESeq2)

# Load data -------------------------------------------------------------------------------------------

load(here::here("results/aberrant_expression/get_gene_count_RSE/gene_counts_rse.rda"))

ref <- dasper:::.ref_load(ref = "/data/references/ensembl/gtf_gff3/v97/Homo_sapiens.GRCh38.97.gtf")

gene_info <- read_delim(here::here("results/aberrant_expression/get_gene_count_RSE/control_1.gene_reads.gct"),
                        delim = "\t", skip = 2) %>% 
  dplyr::select(-control_1)

# Main ------------------------------------------------------------------------------------------------

gene_count_ods <- OutriderDataSet(se = gene_counts_rse)

##### Filter for expressed genes #####

gene_count_ods <- filterExpression(gene_count_ods, 
                                   gtfFile = ref, 
                                   filterGenes = FALSE, 
                                   savefpkm = TRUE)

gene_count_ods <- gene_count_ods[mcols(gene_count_ods)$passedFilter,]

##### Remove the duplicated S2557 sample #####

gene_count_ods <- gene_count_ods[,colData(gene_count_ods)$samp_id_tidy != "S2557"]

##### Run DESeq2 - ATG7 vs rest patients #####

register(MulticoreParam(20))

which_ATG7 <- which(colData(gene_count_ods)[["gene_name"]] == "ATG7")
gene_count_ods$condition <- "control"
gene_count_ods$condition[which_ATG7] <- "case"
gene_count_ods$batch <- gene_count_ods$batch

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
  left_join(gene_info, by = c("gene_id" = "Name")) %>% 
  as_tibble() %>% 
  dplyr::select(gene_symbol = Description, everything())

# rename for save()
ATG7_deseq2_res <- res

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
  left_join(gene_info, by = c("gene_id" = "Name")) %>% 
  as_tibble() %>% 
  dplyr::select(gene_symbol = Description, everything())

COL6A_deseq2_res <- res

# Save data ---------------------------------------------------------------

save(ATG7_deseq2_res, file = here::here("results/aberrant_expression/DESeq2/ATG7_deseq2_res.rda"))

save(COL6A_deseq2_res, file = here::here("results/aberrant_expression/DESeq2/COL6A_deseq2_res.rda"))
