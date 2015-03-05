# Circuitry and Dynamics of Human Transcription Factor Regulatory Networks

Shane Neph
Andrew B. Stergachis
Alex Reynolds
Richard Sandstrom
Elhanan Borenstein
John A. Stamatoyannopoulos

## Highlights

- Extensive transcription factor regulatory networks for 41 human cell and
  tissue types

- Regulatory networks are highly cell selective and expose regulators of
  cellular identity

- Network analysis identifies cell-selective functions for commonly expressed
  regulators

- The circuitry of human transcription factor networks mirrors living neuronal
  networks

## Summary

The combinatorial cross-regulation of hundreds of sequence-specific
transcription factors (TFs) defines a regulatory network that underlies
cellular identity and function. Here we use genome-wide maps of in vivo DNaseI
footprints to assemble an extensive core human regulatory network comprising
connections among 475 sequence-specific TFs and to analyze the dynamics of
these connections across 41 diverse cell and tissue types. We find that human
TF networks are highly cell selective and are driven by cohorts of factors that
include regulators with previously unrecognized roles in control of cellular
identity. Moreover, we identify many widely expressed factors that impact
transcriptional regulatory networks in a cell-selective manner. Strikingly, in
spite of their inherent diversity, all cell-type regulatory networks
independently converge on a common architecture that closely resembles the
topology of living neuronal networks. Together, our results provide an
extensive description of the circuitry, dynamics, and organizing principles of
the human TF regulatory network.

## Files

All network information can be extracted by running: 

  `$ tar -xzf networks.v12032013.tgz`

Each directory (cell type, n=41) includes a file (genes-regulate-genes.txt)
containing all of the TF-to-TF regulatory interactions mapped within that cell
type.

Regulatory interactions were identified in each cell type by scanning the
proximal DHSs (+/-5kb from the canonical transcriptional start site) of each
transcription factor gene DNaseI footprints corresponding to the recognition
sequences of known TFs.

The file genes-regulate-genes.txt contains two columns of gene names. The
gene in the first column is bound (ie contains a DNaseI footprint in proximal
regulatory DNA) corresponding to the second column's TF gene product.
