#' GroupSequences
#'
#'Calculate Hamming distance of a nucleotide sequence alignment and allocate
#' groups according to a threshold.
#'
#' @param aln Nucleotide sequence alignment of class DNAbin or phyDat
#' @param snp_threshold Integer (n>=0) indicating the maximum within-group
#' Hamming distance.
#'
#' @returns A tibble with columns sequence_name and sequence_group.
#' @export
#'
#' @examples
GroupSequences <- function(aln, snp_threshold = 0) {

  if (!inherits(aln, "DNAbin"))
    aln <- ape::as.DNAbin(aln)

  aln <- as.matrix(aln)

  n <- nrow(aln)
  tips <- rownames(aln)

  hd <- hamming_matrix(aln)

  adj <- hd <= snp_threshold
  diag(adj) <- FALSE

  group <- rep(NA_integer_, n)
  gid <- 0L

  for (i in seq_len(n)) {
    if (!is.na(group[i])) next
    gid <- gid + 1L
    stack <- i

    while (length(stack)) {
      v <- stack[[1]]
      stack <- stack[-1]
      if (!is.na(group[v])) next
      group[v] <- gid
      stack <- c(stack, which(adj[v, ]))
    }
  }

  data.frame(
    sequence_name = tips,
    sequence_group = group,
    stringsAsFactors = FALSE
  )
}



################################### Demo #######################################
#library(ape)
#data(woodmouse)
#GroupSequences(woodmouse, snp_threshold = 10)

#################################### END #######################################
################################################################################
