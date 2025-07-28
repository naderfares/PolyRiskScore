# Project Overview

## Purpose
The Polygenic Risk Score Knowledge Base (PRSKB) is a web application and command-line tool for calculating polygenic risk scores using GWAS data from the NHGRI-EBI Catalog.

## Core Components
1. **Web Application**: Node.js/Express server with static HTML/JS frontend
2. **Command-Line Interface**: Bash/Python CLI tool for bulk calculations  
3. **Database Integration**: MySQL database storing GWAS associations and study data

## Key Features
- Calculate polygenic risk scores across multiple traits
- Support for VCF and text input files
- Multiple reference genomes (hg17, hg18, hg19, hg38)
- Population-specific calculations (AFR, AMR, EAS, EUR, SAS)
- Web interface limited to 500 studies, CLI for larger analyses
- JSON and TSV output formats
- Built-in data visualization and comparison tools

## Database Contents (as of March 2022)
- 250,134 variant associations
- 125,433 unique SNPs
- 20,798 study/trait combinations
- 10,366 GWA study identifiers
- 3,463 PubMed identifiers
- Monthly automatic updates from GWAS Catalog