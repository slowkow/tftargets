# tftargets

This package contains two datasets provided by TRED and ITFP:

* `TRED`: Predicted and known human transcription factor targets. (Source:
  https://cb.utdallas.edu/cgi-bin/TRED/)

* `ITFP`: Predicted human transcription factor targets. (Source:
  http://itfp.biosino.org/itfp/)

* `ENCODE`: Putative human transcription factor targets based on ChIP-seq data
  from ENCODE. (Source: http://itfp.biosino.org/itfp/)

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
# Gene symbols used on the TRED website.
> TRED[["STAT3"]]
  [1] "A2M"      "B3GAT3"   "BCL2"     "BCL2L1"   "BIRC5"    "CCL2"     "CCND1"    "CCND3"    "CDKN1A"   "CEBPB"   
 [11] "CSRP1"    "CYP19A1"  "EHHADH"   "FASN"     "FCGR1A"   "FOS"      "HMOX1"    "HSPCA"    "HSPCB"    "IGF1"    
 [21] "IL10"     "IL2RA"    "IL6"      "IL6ST"    "IRF1"     "JAK3"     "JUN"      "KIAA0146" "LBP"      "MCL1"    
 [31] "MIA2"     "MUC1"     "MUC4"     "MYC"      "NOL3"     "NOS2A"    "OSM"      "OXTR"     "PBF"      "PIM1"    
 [41] "PRF1"     "REG1A"    "RORA"     "SEC6L1"   "SOCS1"    "SOCS3"    "SOS1"     "STRA13"   "TIMP1"    "TIMP3"   
 [51] "TLR2"     "TNF"      "TNFRSF5"  "TNFRSF6"  "TNFRSF8"  "TRH"      "VEGF"     "VIP" 

# Gene symbols used on the ITFP website.
> ITFP[["STAT3"]]
 [1] "FIGNL1"   "NCOR1"    "SUV420H1"

# Entrez Gene IDs.
> head(ENCODE[["STAT3"]], 100)
  [1]  23  31  35  40  81  90  93  98 100 104 105 111 114 118 119 135 147 150 159 160 161 174 178 210 224 238 257 259 267
 [30] 272 273 286 287 307 313 320 321 323 328 333 351 368 369 378 402 408 412 419 421 432 444 463 467 472 473 482 491 495
 [59] 529 534 550 571 577 581 586 593 596 597 598 602 622 627 631 636 637 640 651 658 667 669 687 694 695 714 740 752 753
 [88] 770 773 779 780 781 783 788 800 805 811 814 817 821
```
