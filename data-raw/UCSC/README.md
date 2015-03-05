# Transcription Factor ChIP-seq (161 factors) from ENCODE with Factorbook Motifs

<http://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=wgEncodeRegTfbsClusteredV3>

# Description

This track shows regions of transcription factor binding derived from a large
collection of ChIP-seq experiments performed by the ENCODE project, together
with DNA binding motifs identified within these regions by the ENCODE
Factorbook repository.

Transcription factors (TFs) are proteins that bind to DNA and interact with RNA
polymerases to regulate gene expression. Some TFs contain a DNA binding domain
and can bind directly to specific short DNA sequences ('motifs'); others bind
to DNA indirectly through interactions with TFs containing a DNA binding
domain. High-throughput antibody capture and sequencing methods (e.g. chromatin
immunoprecipitation followed by sequencing, or 'ChIP-seq') can be used to
identify regions of TF binding genome-wide. These regions are commonly called
ChIP-seq peaks.

ENCODE TFBS ChIP-seq data were processed using the computational pipeline
developed by the ENCODE Analysis Working Group to generate uniform peaks of TF
binding. Peaks for 161 transcription factors in 91 cell types are combined here
into clusters to produce a summary display showing occupancy regions for each
factor and motif sites within the regions when identified. Additional views of
the underlying ChIP-seq data and documentation on the methods used to generate
it are available from the ENCODE Uniform TFBS track.

# Display Conventions

A gray box encloses each peak cluster of transcription factor occupancy, with
the darkness of the box being proportional to the maximum signal strength
observed in any cell line contributing to the cluster. The HGNC gene name for
the transcription factor is shown to the left of each cluster. Within a
cluster, a green highlight indicates the highest scoring site of a
Factorbook-identified canonical motif for the corresponding factor. (NOTE:
motif highlights are shown only in browser windows of size 50,000 bp or less,
and their display can be suppressed by unchecking the highlight motifs box on
the track configuration page). Arrows on the highlight designate the matching
strand of the motif.

The cell lines where signal was detected for the factor are identified by
single-letter abbreviations shown to the right of the cluster. The darkness of
each letter is proportional to the signal strength observed in the cell line.
Abbreviations starting with capital letters designate ENCODE cell types
identified for intensive study - Tier 1 and Tier 2 - while those starting with
lowercase letters designate Tier 3 cell lines.

Click on a peak cluster to see more information about the TF/cell assays
contributing to the cluster, the cell line abbreviation table, and details
about the highest scoring canonical motif in the cluster.

# Methods

Peaks of transcription factor occupancy from uniform processing of ENCODE
ChIP-seq data by the ENCODE Analysis Working Group were filtered to exclude
datasets that did not pass the integrated quality metric (see "Quality Control"
section of Uniform TFBS) and then were clustered using the UCSC hgBedsToBedExps
tool. Scores were assigned to peaks by multiplying the input signal values by a
normalization factor calculated as the ratio of the maximum score value (1000)
to the signal value at one standard deviation from the mean, with values
exceeding 1000 capped at 1000. This has the effect of distributing scores up to
mean plus one 1 standard deviation across the score range, but assigning all
above to the maximum score. The cluster score is the highest score for any peak
contributing to the cluster.

The Factorbook motif discovery and annotation pipeline uses the MEME-ChIP and
FIMO tools from the MEME software suite in conjunction with machine learning
methods and manual curation to merge discovered motifs with known motifs
reported in Jaspar and TransFac. Motif identifications reported in Wang et al.
2012 (below) were supplemented in this track with more recent data (derived
from newer ENCODE datasets - Jan 2011 through Mar 2012 freezes), provided by
the Factorbook team. Motif identifications from all datasets were merged, with
the most significant value (qvalue) reported being picked when motifs were
duplicated in multiple cell lines. The scores for the selected best-scoring
motif sites were then transformed to -log10.

# Release Notes

Release 4 (February 2014) of this track adds display of the Factorbook motifs.
Release 3 (August 2013) added 124 datasets (690 total, vs. 486 in Release 2),
representing all ENCODE TF ChIP-seq passing quality assessment through the
ENCODE March 2012 data freeze. The peaks used to generate these clusters were
called with less stringent thresholds than used during the January 2011 uniform
processing shown in Release 2 of this track. The contributing datasets are
displayed as individual tracks in the ENCODE Uniform TFBS track, which is
available along with the primary data tracks in the ENC TF Binding Supertrack
page. The clustering for V3/V4 is based on the transcription factor target, and
so differs from V2 where clustering was based on antibody.

For the V3/V4 releases, a new track table format, 'factorSource' was used to
represent the primary clusters table and downloads file,
wgEncodeRegTfbsClusteredV3. This format consists of standard BED5 fields (see
File Formats) followed by an experiment count field (expCount) and finally two
fields containing comma-separated lists. The first list field (expNums)
contains numeric identifiers for experiments, keyed to the
wgEncodeRegTfbsClusteredInputsV3 table, which includes such information as the
experiment's underlying Uniform TFBS table name, factor targeted, antibody
used, cell type, treatment (if any), and laboratory source. The second list
field (expScores) contains the scores for the corresponding experiments. For
convenience, the file downloads directory for this track also contains a BED
file, wgEncodeRegTfbsClusteredWithCellsV3, that lists each cluster with the
cluster score followed by a comma-separated list of cell types.

# Credits

This track shows ChIP-seq data from the Myers Lab at the HudsonAlpha Institute
for Biotechnology and by the labs of Michael Snyder, Mark Gerstein, Sherman
Weissman at Yale University, Peggy Farnham at the University of Southern
California, Kevin Struhl at Harvard, Kevin White at the University of Chicago,
and Vishy Iyer at the University of Texas, Austin. These data were processed
into uniform peak calls by the ENCODE Analysis Working Group pipeline developed
by Anshul Kundaje The clustering of the uniform peaks was performed by UCSC.
The Factorbook motif identifications and localizations (and valuable assistance
with interpretation) were provided by Jie Wang, Bong Hyun Kim and Jiali Zhuang
of the Zlab (Weng Lab) at UMass Medical School.

# References

Gerstein MB, Kundaje A, Hariharan M, Landt SG, Yan KK, Cheng C, Mu XJ, Khurana
E, Rozowsky J, Alexander R et al. Architecture of the human regulatory network
derived from ENCODE data. Nature. 2012 Sep 6;489(7414):91-100. PMID: 22955619

Wang J, Zhuang J, Iyer S, Lin X, Whitfield TW, Greven MC, Pierce BG, Dong X,
Kundaje A, Cheng Y et al. Sequence features and chromatin structure around the
genomic regions bound by 119 human transcription factors. Genome Res. 2012
Sep;22(9):1798-812. PMID: 22955990; PMC: PMC3431495

Wang J, Zhuang J, Iyer S, Lin XY, Greven MC, Kim BH, Moore J, Pierce BG, Dong
X, Virgil D et al. Factorbook.org: a Wiki-based database for transcription
factor-binding data generated by the ENCODE consortium. Nucleic Acids Res. 2013
Jan;41(Database issue):D171-6. PMID: 23203885; PMC: PMC3531197

# Data Release Policy

While primary ENCODE data was subject to a restriction period as described in
the ENCODE data release policy, this restriction does not apply to the
integrative analysis results, and all primary data underlying this track have
passed the restriction date. The data in this track are freely available.
