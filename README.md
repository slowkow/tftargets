# tftargets

This package contains two datasets provided by TRED and ITFP:

* `TRED`: Predicted and known human transcription factor targets. (Source:
https://cb.utdallas.edu/cgi-bin/TRED/)

* `ITFP`: Predicted human transcription factor targets. (Source:
itfp.biosino.org/itfp/)

Install it from github with:
  
```{r}
devtools::install_github("slowikow/tftargets")
```

List the targets of a transcription factor called `STAT3`:

```{r}
> TRED[["STAT3"]]
 [1] "A2M"      "B3GAT3"   "BCL2"     "BCL2L1"   "BIRC5"    "CCL2"     "CCND1"    "CCND3"    "CDKN1A"   "CEBPB"   
[11] "CSRP1"    "CYP19A1"  "EHHADH"   "FASN"     "FCGR1A"   "FOS"      "HMOX1"    "HSPCA"    "HSPCB"    "IGF1"    
[21] "IL10"     "IL2RA"    "IL6"      "IL6ST"    "IRF1"     "JAK3"     "JUN"      "KIAA0146" "LBP"      "MCL1"    
[31] "MIA2"     "MUC1"     "MUC4"     "MYC"      "NOL3"     "NOS2A"    "OSM"      "OXTR"     "PBF"      "PIM1"    
[41] "PRF1"     "REG1A"    "RORA"     "SEC6L1"   "SOCS1"    "SOCS3"    "SOS1"     "STRA13"   "TIMP1"    "TIMP3"   
[51] "TLR2"     "TNF"      "TNFRSF5"  "TNFRSF6"  "TNFRSF8"  "TRH"      "VEGF"     "VIP" 

> ITFP[["STAT3"]]
[1] "FIGNL1"   "NCOR1"    "SUV420H1"
```
