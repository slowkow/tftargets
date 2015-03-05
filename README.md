# tftargets

This package contains the following datasets:

* `TRED`: Predicted and known human transcription factor targets. (Source:
  https://cb.utdallas.edu/cgi-bin/TRED/)

* `ITFP`: Predicted human transcription factor targets. (Source:
  http://itfp.biosino.org/itfp/)

* `ENCODE`: Putative human transcription factor targets based on ChIP-seq data
  from ENCODE. (Source:
  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegTfbsClustered/)

Download and load the RData file:

```{r}
# install.packages("RCurl")
library(RCurl)
download.file(
  url = "https://raw.githubusercontent.com/slowkow/tftargets/master/data/tftargets.RData",
  destfile = "tftargets.RData",
  method = "curl"
)
load("tftargets.RData")
```

List the targets of a transcription factor called `STAT3`:

```{r}
# Entrez Gene IDs.
> TRED[["STAT3"]]
 [1]      2    332    355    595    596    598    896    943    958   1026   1051
[12]   1401   1588   1962   2194   2209   2353   3082   3162   3320   3326   3479
[23]   3559   3572   3586   3659   3718   3725   3929   4170   4582   4585   4609
[34]   4843   5008   5021   5292   5551   5967   6095   6347   6654   7076   7078
[45]   7097   7124   7200   7422   7432   8651   8996   9021  11336  23514  26229
[56]  27151  55893 117153 201254

# Gene symbols used on the ITFP website.
> ITFP[["STAT3"]]
 [1] "FIGNL1"   "NCOR1"    "SUV420H1"

# Entrez Gene IDs.
> head(ENCODE[["STAT3"]], 100)
  [1]  23  31  35  40  81  90  93  98 100 104 105 111 114 118 119 135 147 150 159
 [20] 160 161 174 178 210 224 238 257 259 267 272 273 286 287 307 313 320 321 323
 [39] 328 333 351 368 369 378 402 408 412 419 421 432 444 463 467 472 473 482 491
 [58] 495 529 534 550 571 577 581 586 593 596 597 598 602 622 627 631 636 637 640
 [77] 651 658 667 669 687 694 695 714 740 752 753 770 773 779 780 781 783 788 800
 [96] 805 811 814 817 821
```
