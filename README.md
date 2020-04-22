# Phylogenetic Tree Constructor Documentation

Author: **Ryan Petit**

Date: **21 April 2020**

---

## Development Environment Notes

This program is compatible with the `Python 3.7.6 64-bit` interpreter using the `Anaconda3` distribution. 

It was developed using Microsoft's `Visual Studio Code` IDE within the `macOS Catalina (Version 10.15.3)` operating system.

---

## Dependencies

Please ensure that the following Python modules are installed and up-to-date before process execution:

> numpy: `conda install -c anaconda numpy`

> Bio (Biopython): `conda install -c anaconda biopython`

> matplotlib: `conda install -c anaconda matplotlib`

---

## Process Execution Instructions

This program can be executed using the following terminal command while inside of the `final` directory:

> ```python3 phylogenesis.py``` 

---

## Initialization File Instructions

By default, this program will be configured with the following settings: 

* ALIGNMENT_TYPE: `g`
* SEQUENCE_TYPE: `aa`
* SCORING_MATRIX: `EPAM250`
* MATCH: `5`
* MISMATCH: `-3`
* GAP: `-4`
* SEQUENCE_DIRECTORY: `canis_lupus_aa`
* TREE_NAME: `Unrooted_Phylogram`

These configurations can be edited before the process execution by configuring the `initialization.txt` file also found within the `final` directory. 

Please ensure that the initialization file is formatted in the following way:

> ALIGNMENT_TYPE: `l` (local) or `g` (global)

> SEQUENCE_TYPE: `aa` (amino acid) or `dna` (DNA/mRNA)

> SCORING_MATRIX: `[final name]`

> MATCH: `[integer]`

> MISMATCH: `[integer]`

> GAP: `[integer]`

> SEQUENCE_DIRECTORY: `[directory path from 'final' directory]`

> TREE_NAME: `[string without whitespace]`

Please ensure these values are placed on the right side of the corresponding colon character `':'` for each respective parameter and that the precise name of each parameter goes unedited. Problematic custom configurations will be ignored and will result in default values being used as necessary.

If a parameter is removed from the initialization file, the default value for the parameter will be used during process execution. 

---

## Sequence Batching Instructions

Place all desired FASTA sequences within a single directory. Specify the name of this directory as the value of the `SEQUENCE_DIRECTORY` key inside of the `initialization.txt` file. 

---

## Program Output

The program displays a pop-up MatPlotLib window of the final phylogenetic tree once it has been generated from the multiple sequence alignment.

Each program run generates an output folder within the `ProcessSummaries` directory. The folder is titled with the exact date and time in which the process begins executing. It contains both the `PairwiseAlignments` directory and the `PhyloXML` directory. 

The `PairwiseAlignments` directory is progressively populated with the resultant alignments from the hamming distance matrix population phase of the application. 

The `PhyloXML` directory contains the exported phylogenetic tree in phyloXML format. 