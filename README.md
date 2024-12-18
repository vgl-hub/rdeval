# rdeval

A general purpose, multithreaded read analysis and manipulation tool.

**rdeval** is a single, fast and exhaustive tool for **summary statistics** and simultaneous \*fa\* (fasta, fastq [.gz]) read file **manipulation**. **rdeval** also allows seamless fasta<>fastq[.gz] conversion.

Typical fast\* metrics include:

- number of reads
- total read length
- average read length
- read N50
- shortest read
- longest read
- coverage
- GC content
- base composition

## Installation

Either download one of the releases or `git clone https://github.com/vgl-hub/rdeval.git --recurse-submodules` and `make -j` in the `rdeval` folder.

## Usage

`rdeval input.[fasta|fastq][.gz] [expected genome size]

To check out all options and flags use `rdeval -h`.

You can test some typical usage with the files in the `testFiles` folder, e.g.:

```
rdeval testFiles/random1.fastq 10 // computes summary statistics, including coverage (expected genome size 10bp)
```

## How to cite

If you use **rdeval** in your work, please cite:

Gfastats: conversion, evaluation and manipulation of genome sequences using assembly graphs

Giulio Formenti, Linelle Abueg, Angelo Brajuka, Nadolina Brajuka, Cristo Gallardo, Alice Giani, Olivier Fedrigo, Erich D. Jarvis

doi: https://doi.org/10.1093/bioinformatics/btac460
