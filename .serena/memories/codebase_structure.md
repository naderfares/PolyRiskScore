# Codebase Structure

## Root Directory
- `index.js` - Main Express server entry point (port 3000)
- `package.json` - Node.js dependencies and scripts
- `CLAUDE.md` - Project instructions for Claude Code
- `README.md` - Comprehensive project documentation

## Web Application (`/static/`)
### Frontend Structure
- `*.html` - Static HTML pages (index, calculate_score, studies, visualize, etc.)
- `css/` - Stylesheets
- `js/` - JavaScript modules organized by MVC pattern:
  - `controllers/` - API controllers for database operations
  - `models/` - Database models using Sequelize
  - `routes/` - API route definitions
  - Frontend utility scripts (jQuery, Bootstrap, custom logic)

### CLI Tool (`/static/downloadables/`)
- `runPrsCLI.sh` - Main bash orchestrator script
- `connect_to_server.py` - Downloads GWAS data from server
- `grep_file.py` - Filters input files by relevant SNPs
- `parse_associations.py` - Organizes study/trait data
- `calculate_score.py` - Performs PRS calculations
- `docs/` - CLI documentation
- `README.md` - CLI-specific instructions

## Database Schema
Key tables include:
- `associations` - GWAS variant associations
- `studies` - Study metadata and citations
- `clumps` - Linkage disequilibrium data by population
- `maf` - Minor allele frequencies by cohort

## Other Directories
- `tables/` - Database export/import files
- `update_database_scripts/` - Scripts for data updates
- `.vscode/`, `.idea/` - IDE configurations