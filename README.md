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

Either download the [latest release](https://github.com/vgl-hub/rdeval/releases) or follow the instructions in the [ReadTheDocs documentation](https://rdeval-documentation.readthedocs.io/en/latest/usage.html#installation).

## Usage

`rdeval input.fa*[.gz]|bam|cram|rd [expected genome size]`

To check out all options and flags use `rdeval -h`.

You can test some typical usage with the files in the `testFiles` folder, e.g.:

```
rdeval testFiles/random1.fastq 10 // computes summary statistics, including coverage (expected genome size 10bp)
```

## CiFi support

**rdeval** now supports **SciFi/CiFi single-cell combinatorial indexing** data.  
Using `--cifi-enzyme`, rdeval performs *in-silico* digestion of reads with a specified restriction enzyme motif.  
Using `--cifi-out-combinations`, rdeval outputs all fragment combinations derived from each digested read.

**Basic usage:**

```bash
# Digest reads using enzyme motif
rdeval --cifi-enzyme DpnII scifi_reads.fastq

# Digest and output all combinations
rdeval --cifi-enzyme DpnII --cifi-out-combinations scifi_reads.fastq
```

## Documentation

Additional documentation for rdeval is available in [ReadTheDocs](https://rdeval-documentation.readthedocs.io/).

## How to cite

**Rdeval** is part of the gfastar tool suite. If you use **rdeval** in your work, please cite:

Evaluation of sequencing reads at scale using rdeval.

Formenti G, Koo B, Sollitto M, Balacco J, Brajuka N, Burhans R, Duarte E, Giani AM, McCaffrey K, Medico JA, Myers EW, Smeds P, Nekrutenko A, Jarvis ED. Bioinformatics. 2025 Jul 22:btaf416

doi: https://doi.org/10.1093/bioinformatics/btaf416
PMID: 40694478
