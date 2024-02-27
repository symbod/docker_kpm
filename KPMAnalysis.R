#!/usr/bin/env Rscript

## Script name: KPM.R
##
## Purpose of script: Key Pathway Miner Analysis
##
## Author: Klaudia Adamowicz
## Email: klaudia.adamowicz@uni-hamburg.de
##
## Date Created: 2024-02-25
##
## Copyright (c) Dr. Tanja Laske, 2024

# Load required libraries --------------------------
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
gc()
suppressPackageStartupMessages({
  required_packages <- c("optparse","argparse","data.table", "tibble", "igraph", "tidyverse", "rjson")
  for(package in required_packages){
    if(!require(package,character.only = TRUE)) install.packages(package)
    library(package, character.only = TRUE)
  }
  if(!require("limma",character.only = TRUE, quietly = TRUE)) BiocManager::install("limma")
  library("limma", character.only = TRUE)
  
  # Install KeyPathwayMineR from github and build vignettes
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  if(!require("KeyPathwayMineR",character.only = TRUE)) devtools::install_github("biomedbigdata/keypathwayminer-R", build_vignettes = TRUE)
  # Load and attach KeyPathwayMineR 
  library("KeyPathwayMineR")
})


# Define Methods --------------------------

#' Method performing limma
#'
#' @param set_condition Vector of experimental design specifying the condition(s)
#' @param set_counts Counts (rows = proteins, columns = samples)
#' @param set_comparisons Vector of comparisons
#'
#' @return Results of the limma analysis
#'
#' @export
perform_de <- function(set_condition, set_counts, set_comparisons){
  # create design matrix
  groupsM <- as.factor(set_condition)
  designM <- model.matrix(~0+groupsM) 
  colnames(designM) <-levels(groupsM) 
  fit <- lmFit(set_counts, designM)
  
  # create contrasts
  contr <- makeContrasts(contrasts = set_comparisons, levels = colnames(coef(fit)))
  
  fit2 <- contrasts.fit(fit, contr)
  ebfit <- eBayes(fit2, trend = TRUE)
  return(ebfit)
}

#' Method to extract results of fit object
#'
#' @param set_fit Fit object of the perform_de method
#' @param set_comparisons Vector of comparisons
#' @param out_dir Output directory
#' @param lfc_up Log fold change threshold (upregulated)
#' @param lfc_down Log fold change threshold (downregulated)
#' @param alpha p-value or p.adjust threshold (depending on padj)
#' @param strict TRUE if < and > should be used for threshold comparison, FALSE if <= and >=
#' @param padj TRUE if alpha should be applied to padj, FALSE if alpha applied to p-value
#' @param logFC_thr TRUE if a threshold for logFC should be set, FALSE if not
#'
#' @return Extracted results from the fit object
#'
#' @export
compare_de_expr <- function(set_fit, set_comparisons, out_dir, lfc_up = 1, lfc_down = -1, alpha = 0.05, 
                            strict = FALSE, padj = TRUE, logFC_thr=TRUE, write_output = TRUE){
  de_results <- list()
  for (i in 1:length(set_comparisons)){
    # get results for each comparison
    top.table <- topTable(set_fit, sort.by = "P", number=Inf, coef=c(i))
    gene_reg <- setDT(top.table, keep.rownames = "gene") # save row names as column
    
    # different threshold settings
    if (logFC_thr){
      if (strict){
        if (padj){
          gene_reg$Change <- ifelse(gene_reg$logFC > lfc_up & gene_reg$adj.P.Val < alpha, 1, ifelse(gene_reg$logFC < lfc_down & gene_reg$adj.P.Val < alpha, 1, 0))
        } else {
          gene_reg$Change <- ifelse(gene_reg$logFC > lfc_up & gene_reg$P.Value < alpha, 1, ifelse(gene_reg$logFC < lfc_down & gene_reg$P.Value < alpha, 1, 0))
        }
      } else {
        if (padj){
          gene_reg$Change <- ifelse(gene_reg$logFC >= lfc_up & gene_reg$adj.P.Val <= alpha, 1, ifelse(gene_reg$logFC <= lfc_down & gene_reg$adj.P.Val <= alpha, 1, 0))
        } else {
          gene_reg$Change <- ifelse(gene_reg$logFC >= lfc_up & gene_reg$P.Value <= alpha, 1, ifelse(gene_reg$logFC <= lfc_down & gene_reg$P.Value <= alpha, 1, 0))
        }
      }
    } else {
      if (strict){
        if (padj){
          gene_reg$Change <- ifelse(gene_reg$adj.P.Val < alpha, 1, 0)
        } else {
          gene_reg$Change <- ifelse(gene_reg$P.Value < alpha, 1, 0)
        }
      } else {
        if (padj){
          gene_reg$Change <- ifelse(gene_reg$adj.P.Val <= alpha, 1, 0)
        } else {
          gene_reg$Change <- ifelse(gene_reg$P.Value <= alpha, 1, 0)
        }
      }
    }
    if (write_output){
      write.table(gene_reg, file = paste0(out_dir,"de_data/","DEdata.", set_comparisons[[i]],".txt"),sep=" ",row.names = FALSE) 
    }
    de_results[[set_comparisons[[i]]]] <- gene_reg
  }
  return(de_results)
}

#' Prepare matrix for differential expression analysis
#'
#' @param counts Counts (rows = proteins, columns = samples)
#' @param md Metadata
#' @param condition_name Name of the condition in metadata
#' @param comparisons Vector of comparisons
#' @param out_dir Output directory
#' @param network_proteins Proteins in the network
#' @param lfc_up Log fold change threshold for upregulated proteins
#' @param lfc_down Log fold change threshold for downregulated proteins
#' @param alpha p-value or p.adjust threshold (depending on padj)
#' @param strict TRUE if < and > should be used for threshold comparison, FALSE if <= and >=
#' @param padj TRUE if alpha should be applied to padj, FALSE if alpha applied to p-value
#' @param logFC_thr TRUE if a threshold for logFC should be set, FALSE if not
#' @param plot TRUE to generate plots, FALSE otherwise
#' @param write_output TRUE to write output files, FALSE otherwise
#'
#' @return Prepared matrix for differential expression analysis
#'
#' @export
prepare_matrix <- function(counts, md, condition_name, comparisons, out_dir, network_proteins,
                           lfc_up = args$logFC_up, lfc_down = args$logFC_dow, alpha = args$alpha, 
                           strict = FALSE, padj = args$p_adj, logFC_thr=args$logFC, plot = TRUE, write_output = TRUE){
  if (write_output){
    # create output directories
    dir.create(out_dir, showWarnings = FALSE) #stops warnings if folder already exists
    dir.create(paste0(out_dir,"de_data/"), showWarnings = FALSE) #stops warnings if folder already exists
    
  }
  
  # get condition vector
  condition <- md[, condition_name]
  
  # perform DE analysis pairwise
  pairwise_de <- perform_de(set_condition = condition, set_counts = counts, set_comparisons = comparisons)
  # extract results
  de_results <- compare_de_expr(set_fit = pairwise_de, set_comparisons = comparisons, out_dir = out_dir, lfc_up = lfc_up, 
                                lfc_down = lfc_down, alpha = alpha, strict = strict, padj = padj, logFC_thr = logFC_thr, write_output = write_output)
  mat = data.frame(matrix(ncol=1,nrow=0, dimnames=list(NULL, c("gene"))))
  for (comp in comparisons) {
    cur_df <- de_results[[comp]][,c("gene","Change")]
    colnames(cur_df) <- c("gene", comp)
    mat = merge (mat, cur_df, all=TRUE)
  }
  
  # filter for only expressed rows 
  mat <- mat[rowSums(mat==0) != (ncol(mat)-1), ]
  # explode id column
  mat <- mat %>% mutate(gene=strsplit(gene, ";")) %>% unnest(gene)
  # filter for proteins only in network
  mat <- subset(mat, gene %in% network_proteins)
  return(mat)
}

#' Loop through KPM analysis with varying parameters
#'
#' @param min_k Minimum value of k
#' @param max_k Maximum value of k
#' @param l_min Minimum value of l
#' @param num_pathways Number of pathways to compute
#' @param ind_mat Indicator matrix
#' @param network Network graph
#'
#' @return Result of KPM analysis
#'
#' @export
loop_kpm <- function(min_k = 5, max_k = 10, l_min = 0, num_pathways, ind_mat, network){
  for(k in min_k:max_k){
    # setup kpm
    kpm_options(computed_pathways = num_pathways,k_min = k)
    # run kpm
    result <- kpm(indicator_matrices = ind_mat, graph = network)
    # check found de genes
    found_de <- sum(!result@configurations[[paste0("K-",k,"-L1-",l_min)]]@pathways$`Pathway-1`@nodes$exception)
    if(found_de == nrow(ind_mat)){
      print(paste0("All DE genes found with ", k, " exception nodes."))
      return(result)
    }
  }
  print(paste0("Only ",found_de, "/", nrow(ind_mat), " found."))
  return(result)
}


save_data <- function(out_dir, filename_prefix, df, full) {
  # Construct file paths
  meta_file_path <- file.path(out_dir, paste0(filename_prefix, "_meta.tsv"))
  graphml_file_path <- file.path(out_dir, paste0(filename_prefix, "_network.graphml"))
  
  # Save data frame
  write.table(df, meta_file_path, sep = "\t", row.names = FALSE)
  
  # Create graphml
  g <- graph_from_data_frame(d = full@configurations[[1]]@union_network@edges, 
                             directed = FALSE, 
                             vertices = result_df)
  write_graph(graph = g, file = graphml_file_path, format = "graphml", prefixAttr = FALSE)
}

save_meta_data <- function(df, cond, meta_data, out_dir, filename_prefix) {
  meta_all <- list()
  for (comp in names(df)[-c(1,2,3)]){
    sample_one <- strsplit(comp, split="-")[[1]][1]
    sample_two <- strsplit(comp, split="-")[[1]][2]
    
    # save meta data into json 
    sample_group <- if (cond == "TimeCond") list("Timepoint", "Condition") else if (cond == "TimeCondLoc") list("Timepoint", "Condition", "Location") else list("Condition")
    meta_all[[comp]] = list(
      samples_group_A = meta_data[meta_data[[cond]] %in% c(sample_one)]$Column_name,
      samples_group_B = meta_data[meta_data[[cond]] %in% c(sample_two)]$Column_name,
      group_A = str_split(sample_one,"_")[[1]],
      group_B = str_split(sample_two,"_")[[1]],
      sample_group = sample_group
    )
  }
  # Add links to the entire meta_all list
  meta_all_links <- list(
    dataframe_file = file.path(out_dir, paste0(filename_prefix, "_meta.tsv")),
    graph_file = file.path(out_dir, paste0(filename_prefix, "_network.graphml")))
  
  # Combine metadata with links
  meta_all_combined <- list(meta_data = meta_all, links = meta_all_links)
  # Write to JSON
  write(toJSON(meta_all_combined), file.path(out_dir, paste0(filename_prefix,".txt.json")))
}

## ------------- Prepare Key Pathway Miner Data ------------------

### Parse arguments -----

# set up arguments
parser <- OptionParser()
parser <- add_option(parser, c("-m","--meta_file"), help="Meta data description file")
parser <- add_option(parser, c("-c","--count_file"), help="Preprocessed count file")
parser <- add_option(parser, c("-n","--network_file"), help="Networkfile from chosen organism")
parser <- add_option(parser, c("-o","--out_dir"), help="Directory for output files", default="")

# Adding new options for thresholds with defaults
parser <- add_option(parser, c("--logFC"), help = "Boolean specifying whether to apply a logFC threshold (TRUE) or not (FALSE)", type = "logical", default = TRUE)
parser <- add_option(parser, c("--logFC_up"), help = "Upper log2 fold change threshold (dividing into upregulated)", type = "numeric", default = 1)
parser <- add_option(parser, c("--logFC_down"), help = "Lower log2 fold change threshold (dividing into downregulated)", type = "numeric", default = -1)
parser <- add_option(parser, c("--p_adj"), help = "Boolean specifying whether to apply a threshold on adjusted p-values (TRUE) or on raw p-values (FALSE)", type = "logical", default = TRUE)
parser <- add_option(parser, c("--alpha"), help = "Threshold for adjusted p-values or p-values", type = "numeric", default = 0.05)

# get command line options, if help option encountered print help and exit
args <- parse_args(parser)
# check if mandatory options are set
check_options <- function(tags){
  for (tag in tags){
    if (is.null(args[tag])){
      print_help(parser)
      stop("Missing mandatory option.", call.=FALSE)
    }
  }
}
check_options(c('meta_file','count_file','network_file'))

# save arguments
meta_file_path <- "proteomics-genevention/example_data/plasma/metadata_input.tsv" #args$meta_file
count_file_path <- "proteomics-genevention/example_data/plasma/normalized_counts.tsv"#args$count_file
out_dir <- "proteomics-genevention/example_data/test/kpm/plasma"#args$out_dir
network_file <- "proteomics-genevention/networks/rat_annotated_PPIs_uniprot.sif"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) #stops warnings if folder already exists

#### Load data --------
meta_data <- fread(meta_file_path)
count_data <- fread(count_file_path)
network <- fread(network_file, header=FALSE)

## Prepare data ----

#### Correct data ----
# remove ref
meta_data <- meta_data[meta_data$ID != "ref",]

# convert timepoint column
meta_data[, Timepoint := as.numeric(Timepoint)]

### Rename Columns ---

# rename from file name to sample name
names(count_data) <- plyr::mapvalues(names(count_data), from = meta_data$Column_name, to = meta_data$Sample_name, warn_missing=FALSE) 
meta_data <- subset(meta_data, meta_data$Sample_name %in% names(count_data))

# remove columns that are not in meta_data
columns_to_keep <- c("Protein.IDs", "Gene.Names", meta_data$Sample_name)
existing_columns <- columns_to_keep[columns_to_keep %in% names(count_data)]
count_data <- count_data[, existing_columns, with = FALSE]

#### Encode Columns ----

mappings <- list(
  Condition = list(setNames(paste0("C", seq_along(unique(meta_data$Condition))), unique(meta_data$Condition))),
  Location = if("Location" %in% names(meta_data) && length(unique(meta_data$Location)) > 1) 
    list(setNames(paste0("L", seq_along(unique(meta_data$Location))), unique(meta_data$Location))),
  Timepoint = list(setNames(paste0("T", seq_along(sort(unique(meta_data$Timepoint)))), sort(unique(meta_data$Timepoint))))
)

write(toJSON(mappings), file.path(out_dir, "Metadata_encoding.txt.json")) 

reverse_mappings <- list(
  Condition = if (!is.null(mappings$Condition)) setNames(names(mappings$Condition[[1]]), unlist(mappings$Condition[[1]])) else NULL,
  Location = if (!is.null(mappings$Location)) setNames(names(mappings$Location[[1]]), unlist(mappings$Location[[1]])) else NULL,
  Timepoint = if (!is.null(mappings$Timepoint)) setNames(names(mappings$Timepoint[[1]]), unlist(mappings$Timepoint[[1]])) else NULL
)

write(toJSON(reverse_mappings), file.path(out_dir, "Metadata_reverse_encoding.txt.json")) 

# rename condition
meta_data[, Condition := mappings$Condition[[1]][Condition]]
# rename timepoint
meta_data[, Timepoint := as.character(Timepoint)]
meta_data[, Timepoint := mappings$Timepoint[[1]][Timepoint]]
# rename location
if(!is.null(mappings$Location)) {
  meta_data[, Location := mappings$Location[[1]][Location]]
}

### Create indicator matrix ----
#### Prepare input data ----

# create count matrix 
counts <- as.data.frame(count_data[,c(meta_data$Sample_name), with=F], )
rownames(counts) <- count_data$Protein.IDs
counts <- counts[rowSums(is.na(counts)) != ncol(counts), ] # remove where full row is NA


# filter network
count_data_proteins <- unique(unlist(strsplit(count_data[[1]], ";")))
network <- network[network[[1]] %in% count_data_proteins, ]
network <- network[network[[3]] %in% count_data_proteins, ]
write.table(network, file.path(out_dir, "filtered_network.sif"), sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
filtered_network_file <- file.path(out_dir, "filtered_network.sif")

network_proteins <- unique(c(network$V1, network$V3))

#### Assign active and inactive cases ----
##### Conditions ----
###### Save comparisons ----
comparisons <- c()
elements <- sort(unique(meta_data$Condition))
for (index_a in 1:(length(elements)-1)){
  for (index_b in (index_a+1):length(elements)){
    comparisons <- c(comparisons, paste0(elements[index_a], "-", elements[index_b]))
  }
}
ind_mat_cond <- prepare_matrix(counts=counts, md=as.data.frame(meta_data), condition_name="Condition", comparisons=comparisons, 
                               network_proteins=network_proteins, out_dir = out_dir, write_output = FALSE)
write.table(ind_mat_cond, file.path(out_dir, "indicator_matrix_cond.tsv"), sep="\t", row.names=FALSE)


##### Conditions and Timepoints ----
# save sample groups
target_column <- ifelse("Location" %in% names(meta_data) && length(unique(meta_data$Location)) > 1, 
                        "TimeCondLoc", 
                        "TimeCond")
# Create the target column based on the presence of 'Location'
meta_data[, (target_column) := if(target_column == "TimeCondLoc") {
  paste(meta_data$Timepoint, meta_data$Condition, meta_data$Location, sep = "_")} else {
    paste(meta_data$Timepoint, meta_data$Condition, sep = "_")}]

###### Save comparisons ----
comparisons <- c()
elements <- sort(unique(meta_data[[target_column]]))
for (index_a in 1:(length(elements)-1)){
  for (index_b in (index_a+1):length(elements)){
    # only if time point or condition are the same!
    split_a = str_split(elements[index_a],"_")[[1]]
    split_b = str_split(elements[index_b],"_")[[1]]
    # Count the number of differences
    differences <- sum(split_a != split_b)
    # Only add to comparisons if exactly one of timepoint, condition, or location is different
    if (differences == 1){
      comparisons <- c(comparisons, paste0(elements[index_a], "-", elements[index_b]))
    }
  }
}

ind_mat_timecond <- prepare_matrix(counts=counts, md=as.data.frame(meta_data), condition_name=target_column, comparisons=comparisons, 
                                   network_proteins=network_proteins, out_dir = out_dir, write_output = FALSE)

write.table(ind_mat_timecond, file.path(out_dir,"indicator_matrix_timecond.tsv"), sep="\t", row.names=FALSE)

#### Separate for continuous and static ----

# Function to check if comparison is within the same condition (and optionally location)
is_within_condition <- function(comp, condition, location = NULL) {
  comp_elements <- str_split(comp, "-")[[1]]
  cond_match <- all(str_detect(comp_elements, paste0("_", condition)))
  time_different <- str_split(comp_elements, "_")[[1]][1] != str_split(comp_elements, "_")[[2]][1]
  if (!is.null(location)) {
    loc_match <- all(str_detect(comp_elements, paste0("_", location)))
    return(cond_match && time_different && loc_match)
  } else {
    return(cond_match && time_different)
  }
}

# Function to check if comparison is across conditions at the same time point (and optionally location)
is_across_conditions <- function(comp, location = NULL) {
  comp_elements <- str_split(comp, "-")[[1]]
  time_same <- str_split(comp_elements, "_")[[1]][1] == str_split(comp_elements, "_")[[2]][1]
  condition_different <- str_split(comp_elements, "_")[[1]][2] != str_split(comp_elements, "_")[[2]][2]
  if (!is.null(location)) {
    loc_match <- all(str_detect(comp_elements, paste0("_", location)))
    return(time_same && condition_different && loc_match)
  } else {
    return(time_same && condition_different)
  }
}

# Get unique conditions
unique_conditions <- unique(meta_data$Condition)

# Initialize list to store results
ind_mats <- list()

# Conditional execution based on target_column
if (target_column == "TimeCondLoc") {
  # Get unique locations
  unique_locations <- unique(meta_data$Location)
  
  # Iterate over each location and condition
  for (loc in unique_locations) {
    for (cond in unique_conditions) {
      # Time Point vs Time Point Inside a Condition for each Location
      ind_mat_condx <- ind_mat_timecond[, c(TRUE, sapply(colnames(ind_mat_timecond)[-1], is_within_condition, condition = cond, location = loc))]
      ind_mat_condx <- ind_mat_condx[rowSums(ind_mat_condx == 0) != (ncol(ind_mat_condx)-1), ]
      t_cond <- data.frame(colSums(ind_mat_condx[,-1]))
      
      # Store results
      ind_mats[[paste0(loc, "_", cond, "_tp")]] <- list(
        ind_mat = ind_mat_condx,
        sum_cols = t_cond
      )
    }
    
    # Condition vs Condition per Time Point for each Location
    ind_mat_cond_loc <- ind_mat_timecond[, c(TRUE, sapply(colnames(ind_mat_timecond)[-1], is_across_conditions, location = loc))]
    ind_mat_cond_loc <- ind_mat_cond_loc[rowSums(ind_mat_cond_loc == 0) != (ncol(ind_mat_cond_loc)-1), ]
    
    # Store results
    ind_mats[[paste0(loc, "_cond")]] <- list(
      ind_mat = ind_mat_cond_loc
    )
  }
} else {
  # If Location is not a factor, perform operations without considering Location
  for (cond in unique_conditions) {
    # Time Point vs Time Point Inside a Condition (without location consideration)
    ind_mat_condx <- ind_mat_timecond[, c(TRUE, sapply(colnames(ind_mat_timecond)[-1], is_within_condition, condition = cond))]
    ind_mat_condx <- ind_mat_condx[rowSums(ind_mat_condx == 0) != (ncol(ind_mat_condx)-1), ]
    t_cond <- data.frame(colSums(ind_mat_condx[,-1]))
    
    # Store results
    ind_mats[[paste0(cond, "_tp")]] <- list(
      ind_mat = ind_mat_condx,
      sum_cols = t_cond
    )
  }
  
  # Condition vs Condition per Time Point (without location consideration)
  ind_mat_condx <- ind_mat_timecond[, c(TRUE, sapply(colnames(ind_mat_timecond)[-1], is_across_conditions))]
  ind_mat_condx <- ind_mat_condx[rowSums(ind_mat_condx == 0) != (ncol(ind_mat_condx)-1), ]
  
  # Store results
  ind_mats[["cond"]] <- list(
    ind_mat = ind_mat_condx
  )
}

# Prepare Key Pathway Miner ---------------------

## Save parameters for KPM ----
algorithm = "Greedy" # ("Greedy", "ACO", "Optimal") The algorithm that will be used to extract the pathways
strategy = "INES" # ("GLONE", "INES") The strategy that will be used to extract pathways

# Gene exceptions K (only used for INES):
use_range_k = FALSE # (Boolean) If TRUE k_values will be ranged (batch run). If false the user only needs to set k_min.
k_min = 5 # (Integer) Starting value of k range or k value if k is not ranged
k_step = NULL # (Integer) How k should be increased within the range
k_max = NULL # (Integer) The maximum k value, i.e. the upper limit of the range

# The case (L) exceptions for the n-th matrix:
use_range_l = FALSE # (Boolean) If TRUE l_values will be ranged (batch run). If false the user only needs to set l_min
l_min = 0 # (Integer) Starting value of l range or l value if l is not ranged
l_step = NULL # (Integer) How l should be increased within the range
l_max = NULL # (Integer) The maximum l value, i.e. the upper limit of the range

num_pathways <- 1

# decide what networks should be saved
save_union = TRUE
save_each_pathway = TRUE

settings::reset(kpm_options)

## Set the options for the run ----
kpm_options(
  execution = "Local", 
  strategy = strategy, 
  algorithm = algorithm,
  computed_pathways = num_pathways,
  use_range_k = use_range_k,
  k_min = k_min,
  k_step = k_step,
  k_max = k_max,
  use_range_l = use_range_l,
  l_min = l_min,
  l_step = l_step,
  l_max = l_max,
  remove_bens = TRUE)


# Run Key Pathway Miner -----------------------

conditions <- sort(unique(meta_data$Condition))

get_subnetwork_data <- function(kpm_result){
  result_df <- data.frame(kpm_result@configurations[[1]]@union_network@nodes)
  result_df["size"] <- sapply(strsplit(as.character(result_df[["group"]]), ", "),length)/20
  result_df["group"] <- NULL
  return(result_df)
}

check_submodules <- function(kpm_result, result_df, ind_mat){
  # create constrained subnetwork
  constrained_network <- data.frame(kpm_result@configurations[[1]]@union_network@edges)
  if (nrow(constrained_network) > 0) {
    constrained_network["placeholder"] <- 1
    write.table(constrained_network[c("source","placeholder","target")], file.path(out_dir,"constrained.sif"), sep="\t", quote = FALSE, row.names=FALSE, col.names=FALSE)
  }
  # check each separate module
  for (comp in colnames(ind_mat)[c(-1)]){
    df <- ind_mat[names(ind_mat) %in% c("gene", comp)]
    df <- df[rowSums(df==0) != (ncol(df)-1), ]
    ### run kpm
    if (nrow(df) > 0 && nrow(constrained_network) > 0){
      result <- loop_kpm(min_k = 5, max_k = 10, l_min = 0, num_pathways = 1,
                         ind_mat=data.frame(df), network=file.path(out_dir,"constrained.sif"))
      result_df[comp] <- result_df$node %in% result@configurations[[1]]@pathways$`Pathway-1`@nodes$node
    } else if (nrow(df) > 0 && !nrow(constrained_network) > 0) {
      result_df[comp] <- result_df$node %in% df[[1]]
    } else {
      result_df[comp] <- FALSE
    }
  }
  return(result_df)
}

kpm_results <- list()

## condition vs condition ----
if (nrow(ind_mat_cond) > 0) {
  df <- data.frame(ind_mat_cond$gene, 1)
  #### Loop until all DE nodes found ###
  result_full <- loop_kpm(min_k = 5, max_k = 10, l_min = 0, num_pathways = 20,
                          ind_mat=data.frame(ind_mat_cond), network=filtered_network_file)
  
  result_df <- get_subnetwork_data(kpm_result = result_full)
  result_df <- check_submodules(kpm_result=result_full, result_df=result_df, ind_mat=ind_mat_cond)
  kpm_results[["cond"]] <- list("cond" = result_df)
  
  save_data(out_dir = out_dir, filename_prefix = "cond", df = result_df, full = result_full)
  save_meta_data(df = result_df, cond = "Condition", meta_data = meta_data, 
                 out_dir = out_dir, filename_prefix = "cond")
} else {
  warning("No data found. Skipping condition vs condition analysis.")
}
## time point condition vs time point condition ----

mean_df = separate_rows(count_data,1,sep = ";")[,1]
for (comp in unique(meta_data[[target_column]])) {
  cols = c(names(count_data)[1], meta_data[meta_data[[target_column]] == comp,][["Label"]])
  sub_df = count_data[, cols, with=FALSE]
  sub_df = separate_rows(sub_df,1,sep = ";")
  mean_df[comp] =rowMeans(sub_df[,-1], na.rm = TRUE)
}
write.table(mean_df, file.path(out_dir,"timepoint_mean.tsv"), sep="\t", row.names=FALSE)

kpm_results[["tp"]] <- list()

for ( case in names(ind_mats) ){
  if (nrow(ind_mats[[case]]$ind_mat) > 0) {
    df <- data.frame(ind_mats[[case]]$ind_mat$gene, 1)
    ### run kpm
    result_full <- loop_kpm(min_k = 5, max_k = 10, l_min = 0, num_pathways = 20,
                            ind_mat=df, network=filtered_network_file)
    
    result_df <- get_subnetwork_data(kpm_result=result_full)
    result_df <- check_submodules(kpm_result=result_full, result_df=result_df, ind_mat=ind_mats[[case]]$ind_mat)
    kpm_results[["tp"]][[case]] <- result_df
    
    ### save data
    save_data(out_dir = out_dir, filename_prefix = paste0("timepoint_",case), 
              df = result_df, full = result_full)
    save_meta_data(df = result_df, cond = target_column, meta_data = meta_data, 
                   out_dir = out_dir, filename_prefix = paste0("timepoint_",case))
  } else {
    warning(paste0("No data found. Skipping ", case, " analysis."))
  }
}

# Overview statistics ----

## Top comparisons ----

calculate_exceptions_and_nodes <- function(data) {
  # Extract only the comparison columns
  comparison_cols <- names(data)[!names(data) %in% c("node", "exception", "size")]
  results <- data.frame(comparison = comparison_cols, exceptions = NA, total_nodes = NA)
  
  for (i in seq_along(comparison_cols)) {
    comparison <- comparison_cols[i]
    # Count exceptions and total nodes for each comparison
    results$exceptions[i] <- sum(data[[comparison]] & data$exception, na.rm = TRUE)
    results$total_nodes[i] <- sum(data[[comparison]], na.rm = TRUE)
  }
  
  return(results)
}

## Initialize an empty dataframe for the final results
final_results <- data.frame(source = character(), comparison = character(), exceptions = integer(), total_nodes = integer())

# Loop through each part of kpm_results
for (source_name in names(kpm_results)) {
  for (sub_source_name in names(kpm_results[[source_name]])) {
    # Apply the function to each part
    temp_results <- calculate_exceptions_and_nodes(kpm_results[[source_name]][[sub_source_name]])
    temp_results$source <- paste(source_name, sub_source_name, sep = "_")
    
    # Combine with the final results
    final_results <- rbind(final_results, temp_results)
  }
}

final_results_sorted <- final_results[order(final_results$exceptions, -final_results$total_nodes), ]

write.table(final_results_sorted, file.path(out_dir,"top-comparisons.tsv"), sep="\t", row.names=FALSE)

## Top IDs ----

count_id_occurrences <- function(kpm_results) {
  # Initialize an empty data frame to store the results
  id_results <- data.frame(ID = character(), ExceptionOccurrence = integer(), TotalOccurrence = integer())
  
  # Loop through each main part and sub-part of kpm_results
  for (main_part in names(kpm_results)) {
    for (sub_part in names(kpm_results[[main_part]])) {
      # Current dataframe
      current_df <- kpm_results[[main_part]][[sub_part]]
      
      # Loop through each row of the current dataframe
      for (i in 1:nrow(current_df)) {
        current_id <- current_df$node[i]
        is_exception <- current_df$exception[i]
        
        # If the ID is not in id_results, add it
        if (!current_id %in% id_results$ID) {
          id_results <- rbind(id_results, data.frame(ID = current_id, ExceptionOccurrence = 0, TotalOccurrence = 0))
        }
        
        # Update TotalOccurrence and ExceptionOccurrence
        id_results$TotalOccurrence[id_results$ID == current_id] <- id_results$TotalOccurrence[id_results$ID == current_id] + 1
        if (is_exception) {
          id_results$ExceptionOccurrence[id_results$ID == current_id] <- id_results$ExceptionOccurrence[id_results$ID == current_id] + 1
        }
      }
    }
  }
  
  return(id_results)
}

# Apply the function to kpm_results
id_occurrences <- count_id_occurrences(kpm_results)

id_occurrences_sorted <- id_occurrences[order(-id_occurrences$ExceptionOccurrence, -id_occurrences$TotalOccurrence), ]

write.table(id_occurrences_sorted, file.path(out_dir,"top-ids.tsv"), sep="\t", row.names=FALSE)
