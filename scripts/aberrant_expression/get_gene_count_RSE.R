library(tidyverse)
library(foreach)
library(doParallel)

# Load data ---------------------------------------------------------------

load("/home/dzhang/projects/RNA_seq_diag_mito/results/tidy_samp_metadata/mito_samp_metadata_tidy.rda")

ref <- dasper:::.ref_load(ref = "/data/references/ensembl/gtf_gff3/v97/Homo_sapiens.GRCh38.97.gtf")

# Functions ---------------------------------------------------------------

add_pt_sex <- function(gene_counts_rse, ref){
  
  sex_genes <- GenomicFeatures::genes(ref, filter = list(gene_id = c("ENSG00000067048", "ENSG00000229807")))
  
  gene_counts_rse$sex <- as.character(NA)
  
  for(i in seq_len(dim(gene_counts_rse)[2])){
    
    sex_genes_cov <- rtracklayer::import(gene_counts_rse$bw_path[i], which = sex_genes, as = "NumericList") %>% sum()
    
    stopifnot(sum(sex_genes_cov) > 0)
    
    gene_counts_rse$sex[i] <- 
      ifelse(sex_genes_cov["ENSG00000229807"] > sex_genes_cov["ENSG00000067048"], 
             "female", "male") %>% 
      unname()
    
  }
  
  return(gene_counts_rse)
  
}

# Main --------------------------------------------------------------------

##### Generate gene count matrices #####

numCores <- nrow(mito_samp_metadata_tidy)
registerDoParallel(numCores)

foreach(i=1:nrow(mito_samp_metadata_tidy)) %dopar% {
  
  samp_id <- mito_samp_metadata_tidy[["samp_id_tidy"]][i]
  samp_bam <- mito_samp_metadata_tidy[["bam_path"]][i]
  
  print(stringr::str_c(Sys.time(), " - ", i, " - ", samp_id))
  
  system(
    stringr::str_c("/tools/RNA-SeQC/rnaseqc.v2.3.4.linux", 
                   " /data/references/ensembl/gtf_gff3/v97/Homo_sapiens.GRCh38.97.genes.gtf", 
                   " ", samp_bam, 
                   " /home/dzhang/projects/ATG7_rob_t_analysis/results/get_gene_count_RSE/",
                   " -s ", samp_id, 
                   " -v -v")
    )
  
}

##### Check matrix is as expected #####

gene_count_test <- read_delim("results/get_gene_count_RSE/control_1.gene_reads.gct", 
                              delim = "\t", skip = 2)

summary(gene_count_test$control_1)

# 22746 genes with a TPM above 0 - looks about right

stopifnot(sum(gene_count_test$control_1 > 0) > 1)
sum(gene_count_test$control_1 > 0)

rm(gene_count_test)

##### Merge all patient gene counts #####

gene_count_paths <- list.files("results/get_gene_count_RSE/", pattern = "gene_reads", full.names = TRUE)

for(i in seq_along(gene_count_paths)){
  
  gene_counts <- read_delim(gene_count_paths[i], 
                            delim = "\t", skip = 2)
  
  if(i == 1){
    
    gene_counts_all <- gene_counts
    
  }else{
    
    gene_counts_all <- gene_counts_all %>% 
      left_join(gene_counts)
    
  }
}

##### Create RSE #####

# get gene info and order by genes 
gene_info <- GenomicFeatures::genes(ref, filter = list(gene_id = gene_counts_all[["Name"]]))

stopifnot(all(names(gene_info) %in% gene_counts_all[["Name"]]))
stopifnot(length(gene_info) %in% nrow(gene_counts_all))

gene_counts_all <- gene_counts_all %>% 
  mutate(Name = Name %>% factor(names(gene_info))) %>% 
  dplyr::arrange(Name)

stopifnot(identical(names(gene_info), gene_counts_all[["Name"]] %>% as.character()))

# get sample info and order by samples
gene_counts_mat <- 
  gene_counts_all %>% dplyr::select(-Name, -Description) %>% 
  as.matrix()

mito_samp_metadata_tidy <- mito_samp_metadata_tidy %>% 
  as_tibble() %>% 
  mutate(samp_id_tidy = samp_id_tidy %>% factor(colnames(gene_counts_mat))) %>% 
  arrange(samp_id_tidy)

stopifnot(identical(mito_samp_metadata_tidy[["samp_id_tidy"]] %>% as.character(), 
                    colnames(gene_counts_mat)))

# convert gene counts into an RSE
gene_counts_rse <- SummarizedExperiment::SummarizedExperiment(rowRanges = gene_info, 
                                                           colData = mito_samp_metadata_tidy,
                                                           assays = list(count = gene_counts_mat))

##### add patient sexes #####

gene_counts_rse <- add_pt_sex(gene_counts_rse, ref)

##### add batch #####

gene_counts_rse$batch <- ifelse(str_detect(gene_counts_rse$bw_path, "mito_add_pos_ctrls"), 2, 1)

# Save data ---------------------------------------------------------------

save(gene_counts_rse, file = "results/get_gene_count_RSE/gene_counts_rse.rda")
