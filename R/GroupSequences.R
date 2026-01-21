#' GroupSequences
#'
#'Calculate Hamming distance of a nucleotide sequence alignment and allocate
#' groups according to a threshold.
#'
#' @param aln Nucleotide sequence alignment of class DNAbin or phyDat
#' @param snp_threshold Integer (n>=0) indicating the maximum within-group
#' Hamming distance.
#'
#' @returns A tibble with columns sequence_name and sequence_group
#' @export
#'
#' @examples
GroupSequences <- function(aln, snp_threshold = 0) {

  if (!inherits(aln, c("DNAbin", "phyDat")))
    stop("aln must be DNAbin or phyDat")

  if (snp_threshold < 0)
    stop("SNP threshold must be non-negative")

  # convert formats
  if (inherits(aln, "phyDat")) {
    aln_mat <- as.matrix(as.DNAbin(aln))
    aln_phy <- aln
  } else {
    aln_mat <- as.matrix(aln)
    aln_phy <- as.phyDat(aln)
  }

  n <- nrow(aln_mat)
  tips <- rownames(aln_mat)

  # pairwise SNP distances
  hd <- ape::dist.hamming(aln_phy)
  hd <- as.matrix(hd) * ncol(aln_mat)

  # adjacency matrix
  adj <- hd <= snp_threshold
  diag(adj) <- FALSE

  # ---- connected components (DFS) ----
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
      nbrs <- which(adj[v, ])
      stack <- c(stack, nbrs[group[nbrs] %in% NA])
    }
  }

  data.frame(
    sequence_name  = tips,
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
