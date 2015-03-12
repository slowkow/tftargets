# make-rdata.R
# March 2, 2015
# Kamil Slowikowski

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
  ITFP[[tf]] <- targets
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

# Discard records without an Entrez Gene ID.
dat <- dat[!is.na(dat$entrez), ]

TRED <- list()
for (tf in unique(dat$Factor)) {
  targets <- sort(unique(dat[dat$Factor == tf, "entrez"]))
  TRED[[tf]] <- targets
}

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
  header = FALSE, stringsAsFactors = FALSE, colClasses = "character"
)
peaks <- GRanges(
  seqnames = Rle(peaks$V1),
  ranges = IRanges(
    start = as.numeric(peaks$V2),
    end = as.numeric(peaks$V3)
  ),
  mcols = data.frame(
    TF = peaks$V4,
    Score = as.numeric(peaks$V5)
  )
)

# Scores are scaled such that 1000 is the mean score + 1 standard deviation.
# Also, scores are capped at 1000.
# Requiring a score of 1000 is probably too strict: you can change it.
hits <- findOverlaps(
  genes,
  peaks[peaks$mcols.Score >= 1000, ]
)
dat <- data.frame(
  Gene = genes[queryHits(hits)]$mcols.entrez_id,
  TF = peaks[subjectHits(hits)]$mcols.TF
)

ENCODE <- list()
for (tf in unique(dat$TF)) {
  targets <- sort(unique(dat[dat$TF == tf, "Gene"]))
  ENCODE[[tf]] <- targets
}

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
    Neph2012[[celltype]][[g1]] <- unname(gene2entrez[targets])
  }
#   allgenes <- sort(unique(c(allgenes, dat$V1, dat$V2)))
}

# Save the workspace ----------------------------------------------------------
save(
  list = c("ITFP", "TRED", "ENCODE", "Neph2012"),
  file = "data/tftargets.RData"
)
