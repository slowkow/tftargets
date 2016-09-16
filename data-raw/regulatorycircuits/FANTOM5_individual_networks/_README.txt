FANTOM5 INDIVIDUAL NETWORKS -- regulatorycircuits.org
=====================================================

Transcriptional gene regulatory networks for human cells and tissues based on promoter and
enhancer maps from FANTOM5 described in our paper [1].

   - 394_individual_networks
     The 394 individual cell type and tissue-specific regulatory networks. The mapping of
	 FANTOM5 samples to these networks is given in Supplementary Table 1.

To use these networks with the Magnum app, download first Network_compendium.zip, which
also includes our 32 high-level FANTOM5 networks and other public networks, from:
http://regulatorycircuits.org.

Replace the empty "394_individual_networks" directory in the Network_compendium with the full
directory of the same name included here. Note, it has to be put at precisely this location,
otherwise Magnum will not find it.


FILE FORMAT
===========

The networks are provided as tab-separated text files with three columns, where each line
defines an edge. For directed, weighted networks the format is:
- column 1: the TF gene
- column 2: the target gene
- column 3: the edge weight

Unweighted networks have only 2 columns.
For undirected networks there is no difference between columns 1 and 2.


REFERENCE
=========

[1] Marbach D, Lamparter D, Quon G, Kellis M, Kutalik Z, Bergmann S. Tissue-specific regulatory circuits reveal variable modular perturbations across complex diseases. Submitted.

Contact: daniel.marb...@gmail.com
(fill in the missing letters)


--
Daniel Marbach
September 28, 2015
