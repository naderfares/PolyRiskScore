# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

The Polygenic Risk Score Knowledge Base (PRSKB) is a web application and command-line tool for calculating polygenic risk scores using GWAS data from the NHGRI-EBI Catalog. The project consists of:

1. **Web Application**: Node.js/Express server serving static HTML/JS frontend
2. **Command-Line Interface**: Bash/Python CLI tool for bulk calculations
3. **Database Integration**: MySQL database for storing GWAS associations and study data

## Development Commands

### Server Development
```bash
# Start development server with auto-reload
npm start

# Alternative start (runs nodemon index.js)
node index.js
```

### CLI Tool Testing
```bash
# Run CLI tool (from static/downloadables/)
./runPrsCLI.sh -f inputFile.vcf -o outputFile.tsv -r hg19 -c 0.05 -p EUR

# CLI interactive menu
./runPrsCLI.sh
```

## Architecture

### Web Application Structure
- **index.js**: Main Express server entry point (port 3000)
- **static/**: Frontend assets and CLI downloadables
  - **js/controllers/**: API controllers for database operations
  - **js/models/**: Database models (Sequelize)
  - **js/routes/**: API route definitions
  - **downloadables/**: CLI tool and documentation

### CLI Tool Components
Located in `static/downloadables/`:
- **runPrsCLI.sh**: Main bash script orchestrator
- **connect_to_server.py**: Downloads GWAS data from server
- **grep_file.py**: Filters input files by relevant SNPs
- **parse_associations.py**: Organizes study/trait data
- **calculate_score.py**: Performs PRS calculations

### Database Schema
The application uses MySQL with tables for:
- **associations**: GWAS variant associations
- **studies**: Study metadata and citations
- **clumps**: Linkage disequilibrium data by population
- **maf**: Minor allele frequencies by cohort

### Key Data Processing Pipeline
1. User uploads VCF/TXT file via web interface or CLI
2. System filters variants against study associations
3. LD clumping removes correlated variants
4. PRS calculated using effect sizes and user genotypes
5. Results include percentile rankings against reference populations

## Important Notes

- CLI requires Python3, PyVCF, filelock, requests modules
- Web interface limited to 500 studies per calculation
- Supports reference genomes: hg17, hg18, hg19, hg38
- Population options: AFR, AMR, EAS, EUR, SAS
- No automated testing framework currently implemented