library(tidyverse)
library(OUTRIDER)

# Load data ---------------------------------------------------------------

load("/home/dzhang/projects/ATG7_rob_t_analysis/results/get_gene_count_RSE/gene_counts_rse.rda")

ref <- dasper:::.ref_load(ref = "/data/references/ensembl/gtf_gff3/v97/Homo_sapiens.GRCh38.97.gtf")

# Functions ---------------------------------------------------------------

OUTRIDER_wrapper <- function(rse, out_dir){
  
  ##### Generate an ODS #####
  
  gene_count_ods <- OutriderDataSet(se = rse)
  
  ##### Filter for expressed genes #####
  
  gene_count_ods <- filterExpression(gene_count_ods, 
                                     gtfFile = ref, 
                                     filterGenes = FALSE, 
                                     savefpkm = TRUE)
  
  png(file = stringr::str_c(out_dir, "/FPKM.png"), 
      width = 8, height = 6, units = "in", res = 300)
  
  p <- plotFPKM(gene_count_ods)
  print(p)

  dev.off()
  
  png(file = stringr::str_c(out_dir, "/expressed_genes.png"), 
      width = 8, height = 6, units = "in", res = 300)
  
  p <- plotExpressedGenes(gene_count_ods)
  print(p)
  
  dev.off()
  
  gene_count_ods <- gene_count_ods[mcols(gene_count_ods)$passedFilter,]
  
  ##### Check correlations pre correction #####
  
  png(file = stringr::str_c(out_dir, "/count_cor_pre_correction.png"), 
      width = 8, height = 6, units = "in", res = 300)
  
  gene_count_ods <- plotCountCorHeatmap(gene_count_ods, 
                                        normalized = FALSE,
                                        nRowCluster = 4, 
                                        colGroups = c("batch", "sex"))
  
  dev.off()
  
  png(file = stringr::str_c(out_dir, "/count_gene_sample_pre_correction.png"), 
      width = 8, height = 6, units = "in", res = 300)
  
  gene_count_ods <- plotCountGeneSampleHeatmap(gene_count_ods, 
                                               normalized = FALSE, 
                                               nRowCluster = 4, 
                                               colGroups = c("batch", "sex"))
  
  dev.off()

  ##### Calculate size factors #####

  gene_count_ods <- estimateSizeFactors(gene_count_ods)

  ##### Mask related samples when fitting the autoencoder #####

  mutation_hgvs_concat <- colData(gene_count_ods)[["mutation_hgvs"]] %>%
    lapply(FUN = str_c, collapse = ";") %>% unlist()

  sampleExclusionMask(gene_count_ods) <- (!is.na(mutation_hgvs_concat) & duplicated(mutation_hgvs_concat))

  ##### Find optimal encoding dimension #####

  # This only needs to be run once
  gene_count_ods <- findEncodingDim(gene_count_ods)
  
  png(file = stringr::str_c(out_dir, "/encoding_dim.png"), 
      width = 8, height = 6, units = "in", res = 300)

  p <- plotEncDimSearch(gene_count_ods)
  print(p)
  
  dev.off()

  ##### Correct using autoencoder ######

  which_ATG7 <- which(colData(gene_count_ods)[["gene_name"]] == "ATG7")

  ATG7_res_all <- tibble()

  for(i in which_ATG7){

    samp_interest_id <- colData(gene_count_ods)[["samp_id_tidy"]][i]
    samp_interest <- i
    samp_rest <- which_ATG7[which_ATG7 != i]

    print(stringr::str_c(Sys.time(), " - ", i, " - ", samp_interest_id))

    # remove other ATG7 samples
    gene_count_to_test <- gene_count_ods[,-samp_rest]

    # fit autoencoder
    gene_count_to_test <- controlForConfounders(gene_count_to_test,
                                                q = metadata(gene_count_ods)[["optimalEncDim"]])

    gene_count_to_test <- computePvalues(gene_count_to_test, alternative="two.sided", method="BY")
    gene_count_to_test <- computeZscores(gene_count_to_test)

    png(file = stringr::str_c(out_dir, "/count_cor_post_correction_", samp_interest_id, ".png"),
        width = 8, height = 6, units = "in", res = 300)

    gene_count_to_test <- plotCountCorHeatmap(gene_count_to_test,
                                              normalized = FALSE ,
                                              nRowCluster = 4,
                                              colGroups = c("batch", "sex"))

    dev.off()

    png(file = stringr::str_c(out_dir, "/count_gene_sample_post_correction_", samp_interest_id, ".png"),
        width = 8, height = 6, units = "in", res = 300)

    gene_count_to_test <- plotCountGeneSampleHeatmap(gene_count_to_test,
                                                 normalized = FALSE,
                                                 nRowCluster = 4,
                                                 colGroups = c("batch", "sex"))

    dev.off()

    # volcano <- plotVolcano(gene_count_to_test, samp_interest_id, basePlot=TRUE)
    # 
    # ggsave(filename = stringr::str_c("volcano_", samp_interest_id, ".png"),
    #        path = out_dir,
    #        plot = volcano,
    #        dpi = 300,
    #        width = 8,
    #        height = 6)

    res <- results(gene_count_to_test, all = TRUE)

    res <- res %>% mutate(samp_of_interest = sampleID == samp_interest_id)

    ATG7_res_all <- ATG7_res_all %>%
      bind_rows(res)

  }

  # add gene symbol to results
  gene_count_test <- read_delim("/home/dzhang/projects/ATG7_rob_t_analysis/results/get_gene_count_RSE/control_1.gene_reads.gct",
                                delim = "\t", skip = 2)

  ATG7_res_all <- ATG7_res_all %>%
    left_join(gene_count_test, by = c("geneID" = "Name")) %>%
    dplyr::select(-control_1)

  ##### Save data #####

  save(ATG7_res_all, file = stringr::str_c(out_dir, "/ATG7_res_all.rda"))
  save(gene_count_ods, file = stringr::str_c(out_dir, "/gene_count_ods.rda"))
  save(gene_count_to_test, file = stringr::str_c(out_dir, "/gene_count_corrected_example.rda"))

  return(ATG7_res_all)

}

# Main --------------------------------------------------------------------

OUTRIDER_wrapper(rse = gene_counts_rse, 
                 out_dir = "/home/dzhang/projects/ATG7_rob_t_analysis/results/OUTRIDER/50_controls/")

# filter for only the second batch of ~20 samples
gene_counts_rse <- gene_counts_rse[, gene_counts_rse$batch == 2]

OUTRIDER_wrapper(rse = gene_counts_rse, 
                 out_dir = "/home/dzhang/projects/ATG7_rob_t_analysis/results/OUTRIDER/20_controls/")
