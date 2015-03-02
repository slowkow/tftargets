# make-rdata.R
# March 2, 2015
# Kamil Slowikowski

# ITFP ------------------------------------------------------------------------
dat <- read.delim(
  "inst/extdata/ITFP/regulated_genes.tsv.gz", stringsAsFactors = FALSE)
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
dat <- read.delim("inst/extdata/TRED/targets.tsv.gz", stringsAsFactors = FALSE)
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

TRED <- list()
for (tf in unique(dat$Factor)) {
  targets <- sort(unique(dat[dat$Factor == tf, "Gene"]))
  TRED[[tf]] <- targets
}

sum(names(TRED) %in% names(ITFP)) / length(TRED) # 0.47
sum(names(ITFP) %in% names(TRED)) / length(ITFP) # 0.03

# Check concordance of TRED and ITFP.
sapply(names(TRED)[names(TRED) %in% names(ITFP)], function(tf) {
  t1 <- TRED[[tf]]
  t2 <- ITFP[[tf]]
  total <- sort(unique(c(t1, t2)))
  both <- sum(total %in% t1 & total %in% t2)
  both / length(total)
})

save(
  list = c("ITFP", "TRED"),
  file = "data/tftargets.RData"
)
