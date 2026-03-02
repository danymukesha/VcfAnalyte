# VcfAnalyte - Genetic Variant Analysis Tool

A C-based tool for parsing, annotating, filtering, and summarizing genetic variants from VCF files.

## Features

- VCF 4.2+ compliant parsing
- SNP, indel, CNV, SV detection
- Functional impact prediction
- Filtering by quality, type, frequency, chromosome
- Summary statistics (Ti/Tv ratio, per-chromosome counts)
- **Novel complexity scoring** for identifying challenging variants
- Tab-delimited and JSON output

## Compilation

```bash
gcc -std=c99 -Wall -Wextra -O2 -o vcfanalyte vcfanalyte.c -lm
```

Or use the Makefile: `make`

## Usage

```bash
# Basic with statistics
./vcfanalyte input.vcf -s

# Filter by type and quality
./vcfanalyte input.vcf -t SNP -q 30 -s

# JSON output
./vcfanalyte input.vcf --json > output.json

# Filter by chromosome
./vcfanalyte input.vcf -c 1 -s
```

## Options

- `-s, --stats`         Print summary statistics
- `-j, --json`          Output in JSON format
- `-q, --min-qual`      Minimum quality score
- `-t, --type`          Variant type (SNP, INDEL, SV, CNV)
- `-f, --min-af`        Minimum allele frequency
- `-c, --chrom`          Filter by chromosome
- `-p, --passed-only`   Only PASS variants
- `-o, --output FILE`   Output file

## Input

Accepts VCF files from file path or stdin: `cat input.vcf | ./vcfanalyte -s`
