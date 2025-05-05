<p align="center"><img src="docs/images/gfastar_logo_thumbnail.png" alt="gfastar_logo_thumbnail" width="100" /></p>

# rdeval

A general purpose, multithreaded read analysis and manipulation tool.

**rdeval** is a single, fast and exhaustive tool for **summary statistics** and simultaneous **manipulation** of sequence read files in fa\*[.gz], bam, and cram formats. **rdeval** also allows seamless file conversion between formats.

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

Either download the [latest release](https://github.com/vgl-hub/rdeval/releases/tag/v0.0.7) or follow the instructions in the [ReadTheDocs documentation](https://rdeval-documentation.readthedocs.io/en/latest/usage.html#installation).

## Usage

`rdeval input.fa*[.gz]|bam|cram|rd [expected genome size]`

To check out all options and flags use `rdeval -h`.

You can test some typical usage with the files in the `testFiles` folder, e.g.:

```
rdeval testFiles/random1.fastq 10 // computes summary statistics, including coverage (expected genome size 10bp)
```

Additional documentation for rdeval is available in [ReadTheDocs](https://rdeval-documentation.readthedocs.io/).

## How to cite

**Rdeval** is part of the gfastar tool suite. If you use **rdeval** in your work, please cite:

Gfastats: conversion, evaluation and manipulation of genome sequences using assembly graphs

Giulio Formenti, Linelle Abueg, Angelo Brajuka, Nadolina Brajuka, Cristo Gallardo, Alice Giani, Olivier Fedrigo, Erich D. Jarvis

doi: https://doi.org/10.1093/bioinformatics/btac460
