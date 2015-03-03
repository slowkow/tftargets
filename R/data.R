#' Transcription factor targets from ITFP.
#'
#' Predicted human transcription factor targets scraped from ITFP.
#'
#' @format A list of gene vectors keyed by transcription factor names.
#' \describe{
#'   \item{JUN}{character vector of genes targeted by JUN}
#'   ...
#' }
#' @source \url{itfp.biosino.org/itfp/}
"ITFP"

#' Transcription factor targets from TRED.
#'
#' Predicted and known human transcription factor targets scraped from TRED.
#'
#' @format A list of gene vectors keyed by transcription factor names.
#' \describe{
#'   \item{JUN}{character vector of genes targeted by JUN}
#'   ...
#' }
#' @source \url{https://cb.utdallas.edu/cgi-bin/TRED/}
"TRED"

#' Transcription factor targets estimated from ENCODE data.
#'
#' Human transcription factor targets estimated by selecting a subset of
#' ChIP-seq peaks with a score of 1000.
#'
#' @format A list of gene vectors keyed by transcription factor names.
#' \describe{
#'   \item{JUN}{integer vector of Entrez Gene IDs for genes targeted by JUN}
#'   ...
#' }
#' @source \url{http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegTfbsClustered/}
"ENCODE"
