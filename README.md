# rdeval
Multithreaded read analysis

## Installation

Either download one of the releases or `git clone https://github.com/vgl-hub/rdeval.git --recursive` and `make -j` in `rdeval` folder.

## Usage

`rdeval input.[fasta|fastq|gfa][.gz] [expected genome size]`

To check out all options and flags use `rdeval -h`.

You can test some typical usage example:

```
rdeval testFiles/* 
```

## How to cite

If you use **rdeval** in your work, please cite:

Gfastats: conversion, evaluation and manipulation of genome sequences using assembly graphs

Giulio Formenti, Linelle Abueg, Angelo Brajuka, Nadolina Brajuka, Cristo Gallardo, Alice Giani, Olivier Fedrigo, Erich D. Jarvis

doi: https://doi.org/10.1093/bioinformatics/btac460
