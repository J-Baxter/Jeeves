#' HammingMatrix
#'
#'Calculate Hamming distance of a nucleotide sequence alignment. This is a
#'dependency for GroupSequences.
#'
#' @param aln Nucleotide sequence alignment of class DNAbin or phyDat
#' @param snp_threshold Integer (n>=0) indicating the maximum within-group
#' Hamming distance.
#'
#' @returns A tibble with columns sequence_name and sequence_group.
#' @export
#'
#' @examples
HammingMatrix <- function(aln) {

  n <- nrow(aln)
  L <- ncol(aln)

  mat <- matrix(0L, n, n)

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      mat[i, j] <- mat[j, i] <-
        sum(aln[i, ] != aln[j, ], na.rm = TRUE)
    }
  }

  mat
}
