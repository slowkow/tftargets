# make-rdata.R
# March 2, 2015
# Kamil Slowikowski

library(readr)
library(rjson)
library(RCurl)

# Convert gene names to other IDs.
geneInfo <- function(query) {
  fields <- "symbol,name,taxid,entrezgene"
  url <- sprintf(
    "mygene.info/v2/query?q=%s&entrezonly=1&species=human&fields=%s",
    query, fields
  )
  res <- fromJSON(getURL(url))
  hits <- res$hits
  if (length(hits) > 1) warning(paste("More than 1 hit for", query))
  if (length(hits) < 1) {
    warning(paste("No hits for", query))
    return(NULL)
  }
  res <- do.call(rbind, lapply(hits, function(hit) as.data.frame(hit)))
  res$query <- query
  if (nrow(res[res$X_score > 1, ]) < 1) {
    warning(paste("No good hits for", query))
    return(NULL)
  }
  res[res$X_score > 1, ]
}

# ITFP ------------------------------------------------------------------------
dat <- read.delim(
  "data-raw/ITFP/regulated_genes.tsv.gz", stringsAsFactors = FALSE)
#   TF_name TF_target     score          target_type
# 1    AAAS     ACSF3 0.0431869          normal gene
# 2    AAAS      APRT 0.0458595          normal gene
# 3    AAAS  C12orf52 0.0437940          normal gene
# 4    AAAS    FBXL19 0.0442063 transcription factor
# 5    AAAS     G6PC3 0.0466908          normal gene
# 6    AAMP      NLE1 0.0431282 transcription factor

# Make a list of vectors.
ITFP <- list()
for (tf in unique(dat$TF_name)) {
  targets <- sort(unique(dat[dat$TF_name == tf, "TF_target"]))
  ITFP[[tf]] <- as.character(targets)
}

# TRED ------------------------------------------------------------------------
dat <- read.delim("data-raw/TRED/targets.tsv.gz", stringsAsFactors = FALSE)
gid <- read.delim("data-raw/TRED/gid_entrez.tsv.gz", stringsAsFactors = FALSE)
#  Factor    Gene             Species Map_Location Best_Promoter Best_Promoter_Quality Best_Binding_Quality   gid
# 1   ESR2  GRIN2D human, Homo sapiens  19q13.1-ter         21058 3.1: refseq,predicted                known 14701
# 2   ESR2  XPMC2H human, Homo sapiens       9q34.3         42252              2: known               likely 29252
# 3   ESR2      C3 human, Homo sapiens 19p13.3-13.2         22660              2: known                known 15701
# 4   ESR2 SMARCA1 human, Homo sapiens         Xq25         44093 3.1: refseq,predicted               likely 30486
# 5   ESR2   MKNK2 human, Homo sapiens      19p13.3         22816 3.1: refseq,predicted                known 15811
# 6   ESR2     OXT human, Homo sapiens        20p13         26006              2: known               likely 18040

# Discard non-human.
# table(dat$Species)
dat <- dat[dat$Species == "human, Homo sapiens", ]
dat <- merge(dat, gid, by = "gid")

sum(is.na(dat$entrez)) # 202

dat[is.na(dat$entrez), c("Factor", "Gene", "Map_Location", "gid", "entrez")]

# Manually enter missing Entrez Gene IDs.
# dat[is.na(dat$entrez), c("Factor", "Gene", "gid", "entrez")]
dat[dat$Gene == "PCNA p120", ]$entrez <- 5111
dat[dat$Gene == "DHFR", ]$entrez <- 1719
dat[dat$Gene == "INSIG1", ]$entrez <- 3638
dat[dat$Gene == "SSX8", ]$entrez <- 280659
dat[dat$Gene == "CYP11B2", ]$entrez <- 1585
dat[dat$Gene == "ERVWE1", ]$entrez <- 30816
dat[dat$Gene == "ACDC", ]$entrez <- 9370
dat[dat$Gene == "TRPM2", ]$entrez <- 7226
dat[dat$Gene == "PIGB", ]$entrez <- 9488
dat[dat$Gene == "TU12B1-TY", ]$entrez <- 51559
dat[dat$Gene == "IKBKG", ]$entrez <- 8517
dat[dat$Gene == "SPRR1A", ]$entrez <- 6698

# Fix erroneous Entrez Gene IDs.
dat[dat$Gene == "c-Myc", ]$entrez <- 4609

# Discard records without an Entrez Gene ID.
dat <- dat[!is.na(dat$entrez), ]

TRED <- list()
for (tf in unique(dat$Factor)) {
  targets <- sort(unique(dat[dat$Factor == tf, "entrez"]))
  TRED[[tf]] <- as.character(targets)
}

# Convert TRED TF names to entrez IDs.
d <- do.call(rbind, lapply(names(TRED), geneInfo))
d <- d[!duplicated(d$query), ]

name2entrez <- c(as.character(d$entrezgene), "7161", "1869")
names(name2entrez) <- c(d$query, "TP73L", "E2F-1")

TRED_entrezids <- name2entrez[names(TRED)]
# TRED_symbols <- names(TRED)
# names(TRED_symbols) <- TRED_entrezids

# ENCODE ----------------------------------------------------------------------
library(GenomicRanges)
genes <- read.delim("data-raw/UCSC/knownGene.txt.gz", header = FALSE)
entrez_ids <- read.delim("data-raw/UCSC/knownToLocusLink.txt.gz", header = FALSE)
genes <- merge(genes, entrez_ids, by = "V1")

genes <- GRanges(
  seqnames = Rle(genes$V2.x),
  ranges = IRanges(
    start = genes$V4,
    end = genes$V5
  ),
  strand = Rle(genes$V3),
  mcols = data.frame(
    ucsc_id = genes$V1,
    entrez_id = genes$V2.y
  )
)
# Add 2000 bases upstream of each gene.
genes <- resize(genes, width(genes) + 2000L, fix = "end")
# Add 200 bases downstream of each gene.
genes <- resize(genes, width(genes) + 200L, fix = "start")

peaks <- read.delim(
  "data-raw/UCSC/wgEncodeRegTfbsClusteredV3.bed.gz",
  header = FALSE, stringsAsFactors = FALSE,
  colClasses = c("character", "numeric", "numeric", "character", "integer",
                 "integer", "character", "character")
)
colnames(peaks) <- c(
  "chrom", "start", "end", "name", "score", "expCount", "expNums", "expScores")
peaks <- GRanges(
  seqnames = Rle(peaks$chrom),
  ranges = IRanges(
    start = peaks$start,
    end = peaks$end
  ),
  mcols = data.frame(
    TF = peaks$name,
    Score = peaks$score,
    expCount = peaks$expCount
  )
)

fibroblast_dnase <- read.delim(
  "~/work/fibroblasts/inst/extdata/ENCFF001UUQ.narrowPeak.gz",
  header = FALSE, stringsAsFactors = FALSE,
  colClasses = c("character", "numeric", "numeric", "character", "integer",
                 "character", "numeric", "numeric", "numeric", "numeric")
)
colnames(fibroblast_dnase) <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "signalValue", "pValue", "qValue", "peak")
f_peaks <- GRanges(
  seqnames = Rle(fibroblast_dnase$chrom),
  ranges = IRanges(
    start = fibroblast_dnase$start,
    end = fibroblast_dnase$end
  ),
  mcols = data.frame(
    Score = fibroblast_dnase$score,
    pValue = fibroblast_dnase$pValue
  )
)

# Overlap fibroblast Dnase-seq peaks with ChIP-seq peaks.
# Filters:
#   - Require a p-value < 1e-8.
peaks_subset <- subsetByOverlaps(
  peaks,
  f_peaks[f_peaks$mcols.pValue > 8, ]
)

# Scores are scaled such that 1000 is the mean score + 1 standard deviation.
# Also, scores are capped at 1000.
# Filters:
#   - Require a score of 1000. This is probably too strict: you can change it.
#   - Require that the peak is found in at least 5 experiments.
hits <- findOverlaps(
  genes,
  peaks_subset[
    peaks_subset$mcols.expCount > 10 & peaks_subset$mcols.Score >= 1000, ]
)
dat <- data.frame(
  Gene = genes[queryHits(hits)]$mcols.entrez_id,
  TF = peaks[subjectHits(hits)]$mcols.TF
)

ENCODE <- list()
for (tf in unique(dat$TF)) {
  targets <- sort(unique(dat[dat$TF == tf, "Gene"]))
  ENCODE[[tf]] <- as.character(targets)
}

# Convert ENCODE TF names to entrez IDs.
d <- do.call(rbind, lapply(names(ENCODE), geneInfo))
d <- d[!duplicated(d$query), ]

name2entrez <- c(
  as.character(d$entrezgene),
  "25942", "830649", "10155", "6720", "11128", "55578"
)
names(name2entrez) <- c(
  d$query,
  "SIN3AK20", "GRp20", "KAP1", "SREBP1", "RPC155", "FAM48A"
)

ENCODE_entrezids <- name2entrez[names(ENCODE)]
# ENCODE_symbols <- names(ENCODE)
# names(ENCODE_symbols) <- ENCODE_entrezids

# Neph2012 --------------------------------------------------------------------
celltypes <- basename(list.dirs("data-raw/Neph2012/human_2013-09-16/"))
celltypes <- celltypes[2:length(celltypes)]

# Convert gene symbols to Entrez Gene IDs.
# (I used mygene.info to do the conversion. See "mygene" at Bioconductor.)
dat <- read.delim(
  "data-raw/Neph2012/human_2013-09-16/gene2entrez.tsv",
  stringsAsFactors = FALSE
)
gene2entrez <- dat$entrezgene
names(gene2entrez) <- dat$query

Neph2012 <- list()
# allgenes <- c()
for (celltype in celltypes) {
  dat <- read.delim(file.path(
    "data-raw/Neph2012/human_2013-09-16", celltype, "genes.regulate.genes"
  ), header = FALSE, stringsAsFactors = FALSE)  
  for (g1 in unique(dat$V1)) {
    targets <- sort(unique(dat[dat$V1 == g1, "V2"]))
    targets <- targets[targets %in% names(gene2entrez)]
    targets <- sort(unique(unname(gene2entrez[targets])))
    Neph2012[[celltype]][[g1]] <- as.character(targets)
  }
#   allgenes <- sort(unique(c(allgenes, dat$V1, dat$V2)))
}

# StringDB --------------------------------------------------------------------
# http://string-db.org/api/psi-mi-tab/interactionsList?identifiers=9606.ENSP00000239849&required_score=900

# library(biomaRt)
# mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# results <- getBM(
#   attributes = c("ensembl_gene_id", "ensembl_peptide_id"),
# #   filters = "ensembl_gene_id",
# #   values = "ENSG00000136997",
#   filters = "ensembl_peptide_id",
#   values = "ENSP00000364895",
#   mart = mart
# )
# results

ensembl <- read_tsv(
  file = "data-raw/string-db/ensembl_gene_protein_transcript.txt.gz")
colnames(ensembl) <- c("gene", "protein", "transcript")
ensembl <- ensembl[ensembl$protein != "",]
# gene_to_protein <- with(ensembl, split(protein, gene))
# hist(sapply(gene_to_protein, length))
# table(sapply(gene_to_protein, length))
protein_to_gene <- with(ensembl, split(gene, protein))
protein_to_gene <- protein_to_gene[sapply(protein_to_gene, length) == 1]
protein_to_gene <- unlist(protein_to_gene)
# x <- one_to_one(ensembl, "protein", "gene", exclude = TRUE)

links <- read_delim(
  file = "data-raw/string-db/9606.protein.links.v10.txt.gz", delim = " ")
links <- subset(links, combined_score > 500)

links$protein1 <- substr(links$protein1, 6, 100)
links$protein2 <- substr(links$protein2, 6, 100)

links <- links[links$protein1 %in% names(protein_to_gene),]
links <- links[links$protein2 %in% names(protein_to_gene),]

links$gene1 <- protein_to_gene[links$protein1]
links$gene2 <- protein_to_gene[links$protein2]

stringdb <- split(links$gene1, links$gene2)
stringdb <- lapply(stringdb, function(x) sort(x))

# stringdb2 <- split(links$gene2, links$gene1)
# stringdb2 <- lapply(stringdb2, function(x) sort(x))
# 
# all.equal(stringdb, stringdb2) # TRUE

# length(stringdb) # 2649

# hist(sapply(stringdb, length))

# Save the workspace ----------------------------------------------------------
save(
  list = c(
    "ITFP",
    "TRED", "TRED_entrezids",
    "ENCODE", "ENCODE_entrezids",
    "Neph2012",
    "stringdb"
  ),
  file = "data/tftargets.RData"
#   file = "~/work/fibroblasts/data/tftargets.RData"
)
