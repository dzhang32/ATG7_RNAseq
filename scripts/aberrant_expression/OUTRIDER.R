library(tidyverse)
library(OUTRIDER)

# Load data ---------------------------------------------------------------

load("results/get_gene_count_RSE/gene_counts_rse.rda")

ref <- dasper:::.ref_load(ref = "/data/references/ensembl/gtf_gff3/v97/Homo_sapiens.GRCh38.97.gtf")

# Main --------------------------------------------------------------------

##### Generate an ODS #####

gene_count_ods <- OutriderDataSet(se = gene_counts_rse)

##### Filter by expression #####

gene_count_ods <- filterExpression(gene_count_ods, 
                                   gtfFile = ref, 
                                   filterGenes = FALSE, 
                                   savefpkm = TRUE)

# plotFPKM(gene_count_ods)
# plotExpressedGenes(gene_count_ods)

gene_count_ods <- gene_count_ods[mcols(gene_count_ods)$passedFilter,]

##### Check correlations pre correction #####

gene_count_ods <- plotCountCorHeatmap(gene_count_ods, normalized = FALSE , nRowCluster = 4)

gene_count_ods <- plotCountGeneSampleHeatmap(gene_count_ods, normalized = FALSE, nRowCluster = 4)

##### Calculate size factors #####

gene_count_ods <- estimateSizeFactors(gene_count_ods)

##### Mask related samples when fitting the autoencoder #####

mutation_hgvs_concat <- colData(gene_count_ods)[["mutation_hgvs"]] %>% 
  lapply(FUN = str_c, collapse = ";") %>% unlist()

sampleExclusionMask(gene_count_ods) <- (!is.na(mutation_hgvs_concat) & duplicated(mutation_hgvs_concat))

##### Find optimal encoding dimension #####

# This only needs to be run once 
# Though if rerunning the script and 
gene_count_ods <- findEncodingDim(gene_count_ods)

plotEncDimSearch(gene_count_ods)

##### Correct using autoencoder #####

# gene_count_ods <- controlForConfounders(gene_count_ods, q = 12)
# 
# gene_count_ods <- plotCountCorHeatmap(gene_count_ods, normalized = FALSE , nRowCluster = 4)
# 
# gene_count_ods <- plotCountGeneSampleHeatmap(gene_count_ods, normalized = FALSE, nRowCluster = 4)

##### Obtain results ######

which_ATG7 <- which(colData(gene_count_ods)[["gene_name"]] == "ATG7")

ATG7_res_all <- tibble()

for(i in which_ATG7){
  
  samp_interest_id <- colData(gene_count_ods)[["samp_id_tidy"]][i]
  samp_interest <- i
  samp_rest <- which_ATG7[which_ATG7 != i]
  
  # remove other ATG7 samples
  gene_count_to_test <- gene_count_ods[,-samp_rest]
  gene_count_to_test <- controlForConfounders(gene_count_to_test, q = 12)
  
  gene_count_to_test <- computePvalues(gene_count_to_test, alternative="two.sided", method="BY")
  gene_count_to_test <- computeZscores(gene_count_to_test)
  
  res <- results(gene_count_to_test, all = TRUE)
  
  res <- res %>% mutate(samp_of_interest = sampleID == samp_interest_id) 
  
  ATG7_res_all <- ATG7_res_all %>% 
    bind_rows(res)
  
}

gene_count_test <- read_delim("results/get_gene_count_RSE/control_1.gene_reads.gct", 
                              delim = "\t", skip = 2)

ATG7_res_all <- ATG7_res_all %>% 
  left_join(gene_count_test, by = c("geneID" = "Name")) %>% 
  dplyr::select(-control_1)
  
# Save data ---------------------------------------------------------------

save(ATG7_res_all, file = "results/OUTRIDER/ATG7_res_all.rda")
save(gene_count_to_test, file = "results/OUTRIDER/gene_count_corrected_example.rda")
save(gene_count_ods, file = "results/OUTRIDER/gene_count_ods.rda")
