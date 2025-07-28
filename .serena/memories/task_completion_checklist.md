# Task Completion Checklist

## No Automated Testing
- **Manual testing required** - no automated test framework configured
- Test web interface functionality through browser
- Test CLI tool with sample data files in `static/` directory
- Verify database connections and queries work correctly

## No Linting/Formatting
- **No ESLint or Prettier** configured
- Manual code review for style consistency
- Follow existing code patterns in the codebase
- Ensure consistent indentation and naming

## Development Server Testing
```bash
# Start development server
npm start

# Verify server starts on port 3000
# Test web interface at http://localhost:3000
# Check browser console for JavaScript errors
```

## CLI Tool Testing
```bash
cd static/downloadables/
# Test with sample files
./runPrsCLI.sh -f ../sample.vcf -o test_output.tsv -r hg19 -c 0.05 -p EUR
# Verify output files are generated correctly
```

## Database Verification
- Ensure database connection works (check `static/js/models/database.js`)
- Test API endpoints return expected data
- Verify MySQL queries execute without errors

## File Permissions
- Ensure CLI scripts have execute permissions: `chmod +x runPrsCLI.sh`
- Check Python scripts can be executed
- Verify file upload functionality works

## Documentation Updates
- Update README.md if functionality changes
- Update CLAUDE.md if development patterns change
- Keep CLI documentation in sync with web app changes