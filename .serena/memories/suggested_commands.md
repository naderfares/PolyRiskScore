# Suggested Development Commands

## Server Development
```bash
# Start development server with auto-reload
npm start

# Alternative start (runs nodemon index.js)
node index.js
```

## CLI Tool Testing
```bash
# Navigate to CLI directory
cd static/downloadables/

# Run CLI tool with parameters
./runPrsCLI.sh -f inputFile.vcf -o outputFile.tsv -r hg19 -c 0.05 -p EUR

# Run CLI interactive menu
./runPrsCLI.sh
```

## System Commands (macOS)
```bash
# File operations
ls          # List directory contents
find        # Search for files
grep        # Search text patterns
cd          # Change directory

# Git operations
git status  # Check repository status
git add     # Stage changes
git commit  # Commit changes
git push    # Push to remote
```

## Package Management
```bash
# Install dependencies
npm install

# Update packages
npm update

# Check for security vulnerabilities
npm audit
```

## Database Operations
- Database connection details are in `static/js/models/database.js`
- Credentials stored in `passwords.js` (gitignored)
- Manual database updates via scripts in `update_database_scripts/`