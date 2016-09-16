# tftargets

[Transcription factors (TFs)][TF] activate and repress target genes. This R package
provides easy access to query a particular TF and find its targets in humans.
The data has been collected from multiple different databases.

[TF]: https://www.khanacademy.org/science/biology/gene-regulation/gene-regulation-in-eukaryotes/a/eukaryotic-transcription-factors

<img src="https://github.com/slowkow/tftargets/raw/master/figures/transcription_factor_image.jpg">

Credit: © KENNETH EWARD/BIOGRAFX/PHOTO RESEARCHERS, INC

## Citation

For now, please provide a link to this github repository:

<https://github.com/slowkow/tftargets>


## Download

Download and load the RData file:

```{r}
# Download the file:

# install.packages("RCurl")
library(RCurl)
download.file(
  url = "https://raw.githubusercontent.com/slowkow/tftargets/master/data/tftargets.RData",
  destfile = "tftargets.RData",
  method = "curl"
)

# Load the file:
load("tftargets.RData")

# View the variables stored in the file:
ls()
 [1] "ENCODE"                     "ENCODE_entrezids"           "ITFP"                      
 [4] "Neph2012"                   "regulatory_circuits"        "regulatory_circuits_weight"
 [7] "stringdb"                   "TRED"                       "TRED_entrezids"            
[10] "TRRUST"                     "TRRUST_TYPE"  
```


## Data

This package contains the following datasets:

* [TRED][tred] (2007)
* [ITFP][itfp] (2008)
* [ENCODE][encode] (2012)
* [Neph2012][neph2012] (2012)
* [TRRUST][trrust] (2015)
* [RegulatoryCircuits][regulatorycircuits] (2016)

[tred]: https://github.com/slowkow/tftargets#tred
[itfp]: https://github.com/slowkow/tftargets#itfp
[encode]: https://github.com/slowkow/tftargets#encode
[neph2012]: https://github.com/slowkow/tftargets#neph2012
[regulatorycircuits]: https://github.com/slowkow/tftargets#regulatorycircuits
[trrust]: https://github.com/slowkow/tftargets#trrust

- - -

### TRED

#### Citation

> Jiang, C., Xuan, Z., Zhao, F. & Zhang, M. Q. TRED: a transcriptional
> regulatory element database, new entries and other development. Nucleic
> Acids Res. 35, D137–40 (2007).
> [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/17202159)

#### Source

<https://cb.utdallas.edu/cgi-bin/TRED/tred.cgi?process=home>

#### Description

Predicted and known human transcription factor targets.

Here we find that TRED claims 59 genes are targeted by [STAT3].

[STAT3]: http://www.genecards.org/cgi-bin/carddisp.pl?gene=STAT3

```{r}
# Entrez Gene IDs.
TRED[["STAT3"]]
 [1]      2    332    355    595    596    598    896    943    958   1026   1051
[12]   1401   1588   1962   2194   2209   2353   3082   3162   3320   3326   3479
[23]   3559   3572   3586   3659   3718   3725   3929   4170   4582   4585   4609
[34]   4843   5008   5021   5292   5551   5967   6095   6347   6654   7076   7078
[45]   7097   7124   7200   7422   7432   8651   8996   9021  11336  23514  26229
[56]  27151  55893 117153 201254
```

- - -

### ITFP

#### Citation

> Zheng, G., Tu, K., Yang, Q., Xiong, Y., Wei, C., Xie, L., Zhu, Y. & Li, Y.
> ITFP: an integrated platform of mammalian transcription factors.
> Bioinformatics 24, 2416–2417 (2008).
> [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/18713790)

#### Source

<http://itfp.biosino.org/itfp/>

#### Description

Predicted human transcription factor targets.

```{r}
# Gene symbols used on the ITFP website.
ITFP[["STAT3"]]
 [1] "FIGNL1"   "NCOR1"    "SUV420H1"
```

- - -

### ENCODE

#### Citation

> ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the
> human genome. Nature 489, 57–74 (2012).
> [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/22955616)

#### Source

<http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegTfbsClustered/>

#### Description

Putative human transcription factor targets based on [ChIP-seq] data from
the Encyclopedia of DNA Elements (ENCODE) Project.

[ChIP-seq]: https://en.wikipedia.org/wiki/ChIP-sequencing

```{r}
# Entrez Gene IDs.
head(ENCODE[["STAT3"]], 100)
  [1]  23  31  35  40  81  90  93  98 100 104 105 111 114 118 119 135 147 150 159
 [20] 160 161 174 178 210 224 238 257 259 267 272 273 286 287 307 313 320 321 323
 [39] 328 333 351 368 369 378 402 408 412 419 421 432 444 463 467 472 473 482 491
 [58] 495 529 534 550 571 577 581 586 593 596 597 598 602 622 627 631 636 637 640
 [77] 651 658 667 669 687 694 695 714 740 752 753 770 773 779 780 781 783 788 800
 [96] 805 811 814 817 821
```

- - -

### Neph2012

#### Citation

> Neph, S., Stergachis, A. B., Reynolds, A., Sandstrom, R., Borenstein, E.
> & Stamatoyannopoulos, J. A. Circuitry and dynamics of human transcription
> factor regulatory networks. Cell 150, 1274–1286 (2012).
> [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/22959076)

#### Source

<http://www.regulatorynetworks.org/>

#### Description

Transcription factor targets discovered by DNaseI footprinting and TF
recognition sequences. Targets include only transcription factors and not
other genes.

```{r}
# Entrez Gene IDs.
Neph2012[["AG10803-DS12374"]][["STAT3"]]
 [1]    466   1386    467    468  22809  22926  11016   1385   9586   1390  10664
[12]   1958   1959   1960   1961   2735   2736   2737 148979   2969   8462   9314
[23]   4149   4150   4609   4800   4801   4802   2494   5076   5080   5453   5454
[34]   6667   6668   6670   6671   6774   7020   7021   7022  29842   7490   7494
[45]  51043   7707  10127
```

- - -

### TRRUST

#### Citation

> Han, H., Shim, H., Shin, D., Shim, J. E., Ko, Y., Shin, J., Kim, H., Cho,
> A., Kim, E., Lee, T., Kim, H., Kim, K., Yang, S., Bae, D., Yun, A., Kim, S.,
> Kim, C. Y., Cho, H. J., Kang, B., Shin, S. & Lee, I. TRRUST: a reference
> database of human transcriptional regulatory interactions. Sci. Rep. 5,
> 11432 (2015).
> [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/26066708)

#### Source

<http://www.grnpedia.org/trrust/>

#### Description

> TRRUST is a manually curated database of human transcriptional regulatory
> network.
> 
> Current version of TRRUST contains 8,015 transcriptional regulatory
> relationships between 748 human transcription factors (TFs) and 1,975 non-TF
> genes, derived from 6,175 pubmed articles, which describe small-scale
> experimental studies of transcriptional regulations. To efficiently search
> for regulatory relationships from over 20 million pubmed articles, we used
> sentence-based text mining approach.
> 
> TRRUST database also provide information of mode of regulation (activation
> or repression). Currently 4,861 (60.6%) regulatory relationships are known
> for mode of regulation.

```{r}
head(TRRUST[["STAT3"]], 100)
  [1] "A2M"      "AKAP12"   "AKT1"     "BCL2"     "BCL2"     "BCL2L1"   "BCL2L1"   "BCL6"     "BIRC5"    "BST2"     "CCL11"    "CCL20"   
 [13] "CCND1"    "CCND1"    "CCND2"    "CCND3"    "CD46"     "CDH1"     "CDK4"     "CDKN1A"   "CDKN1B"   "CFB"      "CFLAR"    "CHI3L1"  
 [25] "CISH"     "COPS5"    "CRP"      "CSRP1"    "CTGF"     "CXCL8"    "CYP19A1"  "CYR61"    "DDIT3"    "DNMT1"    "EGFR"     "ESR2"    
 [37] "ETV6"     "F2R"      "FAAH"     "FAS"      "FAS"      "FGF1"     "FGF2"     "FGG"      "FGL1"     "FLT3"     "FOS"      "GAST"    
 [49] "GFAP"     "HAMP"     "HGF"      "HIF1A"    "HMOX1"    "HP"       "HSPA4"    "HSPB1"    "ICAM1"    "IFNAR1"   "IFNG"     "IKBKE"   
 [61] "IL10"     "IL11"     "IL1RN"    "IL2"      "IL21"     "IL2RA"    "IL6"      "IL6"      "IRF1"     "JAK2"     "JAK3"     "JUNB"    
 [73] "KLF11"    "KRT17"    "LCAT"     "LEP"      "LGALS3BP" "LTBP1"    "MCL1"     "MCL1"     "MDC1"     "MICA"     "MMP1"     "MMP14"   
 [85] "MMP2"     "MMP2"     "MMP3"     "MMP7"     "MMP7"     "MMP9"     "MMP9"     "MUC1"     "MUC4"     "MYC"      "MYC"      "NANOG"   
 [97] "NDUFA13"  "NME1"     "NOSTRIN"  "NOX5"  
```

#### Raw Data

```bash
zcat data-raw/TRRUST/trrust_rawdata.txt.gz | head | column -t
AATF  BAK1    Unknown     22983126
AATF  BAX     Repression  22909821
AATF  BBC3    Unknown     22983126
AATF  CDKN1A  Unknown     17157788
AATF  MYC     Activation  20549547
AATF  TP53    Unknown     17157788
ABL1  BAX     Activation  11753601
ABL1  BCL2    Repression  11753601
ABL1  BCL6    Repression  15509806
ABL1  CCND2   Activation  15509806
```

<img src="https://github.com/slowkow/tftargets/raw/master/figures/TRRUST_histogram.png">

- - -

### Regulatory Circuits

#### Citation

> Marbach, D., Lamparter, D., Quon, G., Kellis, M., Kutalik, Z. & Bergmann, S.
> Tissue-specific regulatory circuits reveal variable modular perturbations
> across complex diseases. Nat. Methods 13, 366–370 (2016).
> [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/26950747)

#### Source

<http://regulatorycircuits.org>

#### Description

> We developed a comprehensive resource of close to 400 cell type- and
> tissue-specific gene regulatory networks for human. Our study shows that
> disease-associated genetic variants often perturb regulatory modules in cell
> types or tissues that are highly specific to that disease.

```{r}
head(regulatory_circuits[["STAT3"]], 100)
  [1] "SURF1"    "ZNF230"   "EIF5"     "ATG4C"    "LYSMD4"   "ZWILCH"   "TFB1M"    "SLC12A7"  "DNAL1"    "PPP1R8"   "SEPT9"    "SDCCAG8" 
 [13] "CMTR1"    "GSAP"     "PPIA"     "CLCN6"    "ZFP69"    "ZFP64"    "RNPC3"    "BRPF1"    "ZKSCAN5"  "ZNF410"   "ASF1B"    "PES1"    
 [25] "TMEM41B"  "F2RL1"    "DARS"     "ZNF24"    "RPL4"     "SYF2"     "AGTPBP1"  "NANOS1"   "ZNF140"   "SEC14L1"  "CHAC1"    "CDC42SE2"
 [37] "LIPG"     "PROS1"    "MIIP"     "DENND1A"  "ADAMTSL2" "TBC1D22B" "PHACTR4"  "TNFAIP2"  "SLC35C1"  "ZNF284"   "NCCRP1"   "ZFYVE16" 
 [49] "TBL1XR1"  "UNC45A"   "TIMM50"   "PRRT1"    "RNF215"   "PAF1"     "SPINT1"   "RABL2B"   "DMWD"     "RIN3"     "PAK2"     "NOTCH4"  
 [61] "INPP5F"   "PSMA8"    "MX2"      "TBC1D7"   "CCDC135"  "ATP2B4"   "HLA-DQA2" "IPO8"     "EID2B"    "OGDH"     "ZFYVE21"  "DDB1"    
 [73] "SEC31A"   "SURF6"    "EXD2"     "KIF3A"    "RPUSD3"   "SYMPK"    "ASB13"    "CASC5"    "RLF"      "LIN54"    "TNXB"     "TRABD"   
 [85] "PHTF2"    "COPS4"    "FAM32A"   "PDLIM4"   "CPSF7"    "ZNF720"   "RBFOX2"   "COA4"     "ATP10A"   "MTMR1"    "TNRC6C"   "TMED4"   
 [97] "BUD31"    "GADD45B"  "MTMR3"    "CDC42EP4"
```

#### Raw Data

Columns:

1. Transcription factor.
2. Target gene.
3. Edge weight.

```bash
zcat data-raw/regulatorycircuits/FANTOM5_individual_networks/394_individual_networks/synoviocyte.txt.gz | head
RAX    PPP2R2A  1.79016453E-3
MYCN   RHOA     1.81311653E-2
TFAP2  RRM1     7.13096624E-3
PRDM4  KPNA2    1.61069158E-2
FOXB1  SCARF2   1.78696733E-3
ATF4   NDUFA11  1.53625527E-3
SPIC   C9orf69  8.60099271E-4
FLI1   CENPU    6.72504942E-3
HNF4A  LHFPL2   1.47391413E-2
STAT3  SURF1    3.14614561E-3
```

<img src="https://github.com/slowkow/tftargets/raw/master/figures/regulatory_circuits_histogram.png">

<img src="https://github.com/slowkow/tftargets/raw/master/figures/regulatory_circuits_weights_histogram.png">

<img src="https://github.com/slowkow/tftargets/raw/master/figures/regulatory_circuits_targets_vs_weight.png">

