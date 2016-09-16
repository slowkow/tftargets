#' Transcription factor targets from TRED.
#'
#' Predicted and known human transcription factor targets scraped from TRED.
#'
#' Jiang, C., Xuan, Z., Zhao, F. & Zhang, M. Q. TRED:
#' a transcriptional regulatory element database, new entries and other
#' development. Nucleic Acids Res. 35, D137–40 (2007).
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/17202159}{PubMed})
#'
#' @docType data
#'
#' @format A list of gene vectors keyed by transcription factor names.
#' \describe{
#'   \item{JUN}{integer vector of Entrez Gene IDs for genes targeted by JUN}
#'   ...
#' }
#'
#' @source \url{https://cb.utdallas.edu/cgi-bin/TRED/}
#'
#' @examples
#' length(TRED)
#' TRED[1]
#' \donttest{hist(sapply(TRED, length))}
"TRED"

#' Transcription factor targets from ITFP.
#'
#' Predicted human transcription factor targets scraped from ITFP.
#'
#' Zheng, G., Tu, K., Yang, Q., Xiong, Y., Wei, C., Xie, L., Zhu,
#' Y. & Li, Y. ITFP: an integrated platform of mammalian transcription
#' factors. Bioinformatics 24, 2416–2417 (2008).
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/18713790}{PubMed})
#'
#' @docType data
#'
#' @format A list of gene vectors keyed by transcription factor names.
#' \describe{
#'   \item{JUN}{character vector of genes targeted by JUN}
#'   ...
#' }
#'
#' @source \url{itfp.biosino.org/itfp/}
#'
#' @examples
#' length(ITFP)
#' ITFP[1]
#' \donttest{hist(sapply(ITFP, length))}
"ITFP"

#' Transcription factor targets estimated from ENCODE data.
#'
#' Human transcription factor targets estimated by selecting a subset of
#' ChIP-seq peaks with a score of 1000.
#'
#' The ChIP-seq data used here is taken from sample \code{ENCFF001UUQ}.
#' There are many other samples that you might want to use instead.
#'
#' ENCODE Project Consortium. An integrated encyclopedia of DNA
#' elements in the human genome. Nature 489, 57–74 (2012).
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/22955616}{PubMed})
#'
#' @docType data
#'
#' @format A list of gene vectors keyed by transcription factor names.
#' \describe{
#'   \item{JUN}{integer vector of Entrez Gene IDs for genes targeted by JUN}
#'   ...
#' }
#' @source \url{http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegTfbsClustered/}
#'
#' @examples
#' length(ENCODE)
#' ENCODE[1]
#' \donttest{hist(sapply(ENCODE, length))}
"ENCODE"

#' Transcription factor targets from Neph et al. 2012.
#'
#' Human transcription factor targets obtained from DNaseI footprinting and TF
#' recognition sequences. Targets only include transcription factors.
#'
#' Neph, S., Stergachis, A. B., Reynolds, A., Sandstrom, R.,
#' Borenstein, E. & Stamatoyannopoulos, J. A. Circuitry and dynamics of human
#' transcription factor regulatory networks. Cell 150, 1274–1286 (2012).
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/22959076}{PubMed})
#'
#' @docType data
#'
#' @format A list of lists keyed by cell type names. Each sublist is keyed by
#' transcription factor names and returns a gene vector.
#' \describe{
#'   \item{JUN}{integer vector of Entrez Gene IDs for genes targeted by JUN}
#'   ...
#' }
#'
#' @source \url{http://www.regulatorynetworks.org/}
#'
#' @examples
#' length(Neph2012)
#' Neph2012[1]
#' \donttest{hist(sapply(Neph2012, length))}
"Neph2012"

#' Transcription factor targets from TRRUST.
#'
#' Han, H., Shim, H., Shin, D., Shim, J. E., Ko, Y., Shin, J.,
#' Kim, H., Cho, A., Kim, E., Lee, T., Kim, H., Kim, K., Yang, S., Bae, D.,
#' Yun, A., Kim, S., Kim, C. Y., Cho, H. J., Kang, B., Shin, S. & Lee, I.
#' TRRUST: a reference database of human transcriptional regulatory
#' interactions. Sci. Rep. 5, 11432 (2015).
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/26066708}{PubMed})
#'
#' @docType data
#'
#' @format A list of gene vectors keyed by transcription factor names.
#' \describe{
#'   \item{JUN}{integer vector of Entrez Gene IDs for genes targeted by JUN}
#'   ...
#' }
#' @source \url{http://www.grnpedia.org/trrust/}
#'
#' @examples
#' length(TRRUST)
#' TRRUST[1]
#' \donttest{hist(sapply(TRRUST, length))}
"TRRUST"

#' Transcription factor targets from Marbach et al. 2016.
#'
#' This data is taken from the \code{synoviocyte.txt} file. There are many
#' other cell types available that you might want to use instead.
#'
#' Marbach, D., Lamparter, D., Quon, G., Kellis, M., Kutalik, Z.
#' & Bergmann, S. Tissue-specific regulatory circuits reveal variable modular
#' perturbations across complex diseases. Nat. Methods 13, 366–370 (2016).
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/26950747}{PubMed})
#'
#' @docType data
#'
#' @format A list of gene vectors keyed by transcription factor names.
#' \describe{
#'   \item{JUN}{character vector of genes targeted by JUN}
#'   ...
#' }
#'
#' @source \url{http://www.regulatorycircuits.org/}
#'
#' @examples
#' length(regulatory_circuits)
#' regulatory_circuits[1]
#' \donttest{hist(sapply(regulatory_circuits, length))}
"regulatory_circuits"
