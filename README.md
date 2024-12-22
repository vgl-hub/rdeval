<p style="text-align:center;"><img src="docs/images/gfastar_logo_thumbnail.png" alt="gfastar_logo_thumbnail" width="250"></p>

# rdeval

A general purpose, multithreaded read analysis and manipulation tool.

**rdeval** is a single, fast and exhaustive tool for **summary statistics** and simultaneous **manipulation** of sequence read files in fa\*[.gz], bam, cram, and formats. **rdeval** also allows seamless file conversion conversion between formats.

Typical metrics include:

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

Either download one of the releases or `git clone https://github.com/vgl-hub/rdeval.git --recursive` and `make -j` in the `rdeval` folder.

## Usage

`rdeval input.fa*[.gz]|bam|cram|rd [expected genome size]`

To check out all options and flags use `rdeval -h`.

You can test some typical usage with the files in the `testFiles` folder, e.g.:

```
rdeval testFiles/random1.fastq 10 // computes summary statistics, including coverage (expected genome size 10bp)
```

## How to cite

**Rdeval** is part of the gfastar tool suite. If you use **rdeval** in your work, please cite:

Gfastats: conversion, evaluation and manipulation of genome sequences using assembly graphs

Giulio Formenti, Linelle Abueg, Angelo Brajuka, Nadolina Brajuka, Cristo Gallardo, Alice Giani, Olivier Fedrigo, Erich D. Jarvis

doi: https://doi.org/10.1093/bioinformatics/btac460
