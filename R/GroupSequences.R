################################################################################
## Script Name:        GroupSequences
## Purpose:            Calculate Hamming distance of a nucleotide sequence
##                     alignment and allocate groups according to a threshold.
## Author:             James Baxter
## Date Created:       2025-08-14
################################################################################


################################### MAIN #######################################
# Main analysis or transformation steps
GroupSequences <- function(aln, snp_threshold = 0){
  require(igraph)
  require(phangorn)
  require(tidyverse)
  require(magrittr)

  # Sanity checks
  if(class(aln) != "DNAbin"){
    stop('Aln must be a DNAbin.')
  }

  if(!is.matrix(aln)){
    stop('Aln must be a matrix.')
  }

  if(snp_threshold < 0){
    stop('SNP threshold must be positive.')
  }


  # Ensure alignment is correctly formatted
  if(class(aln) != 'PhyDat'){
    aln_formatted <- as.phyDat(aln)
  }else{
    aln_formatted <- aln
  }

  # Calculate hamming distance
  hd_normalised <- dist.hamming(aln_formatted) %>%
    as.matrix()
  hd_raw <- hd_normalised * ncol(aln)

  # Obtain groups of sequences for which HD < SNP threshold
  if( any(hd_raw[lower.tri(hd_raw, diag = FALSE)] <= snp_threshold)){
    groups <- which(hd_raw <= snp_threshold,
                    arr.ind = TRUE) %>%
      dplyr::as_tibble(rownames = 'tipnames') %>%
      filter(row !=col) %>%
      dplyr::select(-tipnames) %>%

      # Infer network from HDs
      igraph::graph_from_data_frame(.,
                                    directed = F) %>%

      components() %>%
      getElement('membership') %>%
      stack() %>%
      as_tibble() %>%
      mutate(ind = as.numeric(as.character(ind))) %>%
      mutate(tipnames = map_chr(ind, ~ rownames(aln)[.x])) %>%
      dplyr::select(c(tipnames, values)) %>%
      dplyr::distinct() %>%
      dplyr::rename(sequence_group = values)

  }else{
    warning('all sequences above threshold.')
    groups <- tibble(tipnames = rownames(aln)) %>%
      rowid_to_column(., var = 'sequence_group')
  }


  out <- tibble(tipnames = rownames(aln)) %>%
    left_join(groups) %>%
    mutate(sequence_group =
             ifelse(is.na(sequence_group),
                    max(sequence_group, na.rm = T) + row_number() + 1,
                    sequence_group))


  return(out)
}


################################### Demo #######################################
#library(ape)
#data(woodmouse)
#GroupSequences(woodmouse, snp_threshold = 10)

#################################### END #######################################
################################################################################
