# Load packages -----------------------------------------------------------

library(tidyverse)
library(stringr)

# this script requires this exact version of dasper
# it's likely that the future updates to dasper may break some the code used
# so to avoid this you can install using the commit hash below
# devtools::install_github("dzhang32/dasper@6e87cf1f403b06548f12cff9290029f22a54b781")
library(dasper)

# Load data ---------------------------------------------------------------

load("/home/dzhang/projects/RNA_seq_diag_mito/results/tidy_samp_metadata/mito_samp_metadata_tidy.rda")

blacklist_hg38 <- rtracklayer::import("/data/references/ENCODE_blacklist_v2/hg38-blacklist.v2.bed")

ref <- "/data/references/ensembl/gtf_gff3/v97/Homo_sapiens.GRCh38.97.gtf"
suppressWarnings(expr = {
  ref <- dasper:::.ref_load(ref)
})

# Main --------------------------------------------------------------------

##### Load junctions #####

# junction_load() doesn't currently handle DataFrames for metadata
# as these are incompatible with bind_rows()
# convert to data.frame - will change CharacterList cols to lists
pos_control_metadata <- as.data.frame(mito_samp_metadata_tidy)

junctions <- dasper::junction_load(junction_paths = pos_control_metadata[["junctions_path"]],
                                   metadata = pos_control_metadata,
                                   chrs = c(1:22, "X", "Y", "MT"))

# convert lists back into CharacterLists
for(i in seq_along(SummarizedExperiment::colData(junctions))){
  
  if(is.list(SummarizedExperiment::colData(junctions)[[i]])){
    
    SummarizedExperiment::colData(junctions)[[i]] <- SummarizedExperiment::colData(junctions)[[i]] %>%
      IRanges::CharacterList()
    
  }
  
}

##### Set ATG7 patients as cases #####

# remove technical replicate
junctions <- junctions[, SummarizedExperiment::colData(junctions)[["samp_id_tidy"]] != "S2557"]

# add control identifier for all remaining patients (except ATG7)
# need to generate a separate object for 
patient_metadata <- SummarizedExperiment::colData(junctions)
which_ATG7 <- which(patient_metadata[["gene_name"]] == "ATG7")

SummarizedExperiment::colData(junctions)[["case_control"]] <- "control"
SummarizedExperiment::colData(junctions)[["case_control"]][which_ATG7] <- "case"

##### Process junctions #####

# convert seqlevels to match junctions
GenomeInfoDb::seqlevels(blacklist_hg38) <- GenomeInfoDb::seqlevels(blacklist_hg38) %>%
  str_replace("chr", "")
stopifnot(all(GenomeInfoDb::seqlevels(blacklist_hg38) %in% GenomeInfoDb::seqlevels(junctions)))

# saving this intermediary object for plotting
junctions_normed <- junctions %>% 
  junction_filter(
    count_thresh = c("raw" = 5),
    n_samp = c("raw" = 1),
    width_range = c(20, 1000000),
    regions = blacklist_hg38
  ) %>% 
  junction_annot(ref) %>% 
  junction_filter(types = c("unannotated", "ambig_gene")) %>% 
  junction_norm()

junctions <- junctions_normed %>% 
  junction_score()

##### Process coverage #####

cases <- SummarizedExperiment::colData(junctions_normed)[["case_control"]] == "case"
controls <- SummarizedExperiment::colData(junctions_normed)[["case_control"]] == "control"
bw_path_case <- SummarizedExperiment::colData(junctions_normed)[["bw_path"]][cases]
bw_path_control <- SummarizedExperiment::colData(junctions_normed)[["bw_path"]][controls]

stopifnot(length(bw_path_case) == 5)
stopifnot(length(bw_path_control) == 53)

junctions <- coverage_process(
  junctions, 
  ref, 
  coverage_paths_case = bw_path_case, 
  coverage_paths_control = bw_path_control, 
  bp_param = BiocParallel::MulticoreParam(length(bw_path_control))
)



##### Coverage_score/junction_score #####

outlier_scores <- outlier_process(
  junctions, 
  samp_id_col = "samp_id_tidy", 
  random_state = 32L
)

# Save data ---------------------------------------------------------------

save(junctions_normed, 
     file = "/data/junctions_normed_for_plotting.rda", 
     compress = "gzip")

save(junctions, 
     file = "/data/junctions.rda", 
     compress = "gzip")

save(outlier_scores, 
     file = "/data/outlier_scores.rda", 
     compress = "gzip")