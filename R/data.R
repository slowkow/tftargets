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
#'   \item{JUN}{integer vector of Entrez Gene IDs for genes targeted by JUN}
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

#' Transcription factor targets from Neph et al. 2012.
#'
#' Human transcription factor targets obtained from DNaseI footprinting and TF
#' recognition sequences. Targets only include transcription factors.
#'
#' @format A list of lists keyed by cell type names. Each sublist is keyed by
#' transcription factor names and returns a gene vector.
#' \describe{
#'   \item{JUN}{integer vector of Entrez Gene IDs for genes targeted by JUN}
#'   ...
#' }
#' @source \url{http://www.regulatorynetworks.org/}
"Neph2012"

#' Transcription factor targets from Neph et al. 2012.
#'
#' @format A list of gene vectors keyed by transcription factor names.
#' \describe{
#'   \item{JUN}{integer vector of Entrez Gene IDs for genes targeted by JUN}
#'   ...
#' }
#' @source \url{http://www.grnpedia.org/trrust/}
"TRRUST"

#' Transcription factor targets from Marbach et al. 2016.
#'
#' @format A list of gene vectors keyed by transcription factor names.
#' \describe{
#'   \item{JUN}{character vector of genes targeted by JUN}
#'   ...
#' }
#' @source \url{http://www.regulatorycircuits.org/}
"regulatory_circuits"
