# Code Style and Conventions

## JavaScript Style
- **No formal linting** or formatting tools configured
- Uses **ES5/ES6 mixed syntax** throughout codebase
- **jQuery-based** frontend with traditional script organization
- **Callback-style** async patterns (not async/await in most places)
- **CommonJS** require/module.exports for Node.js modules

## File Naming
- **camelCase** for JavaScript files and functions
- **kebab-case** for HTML files
- **PascalCase** for model/controller constructors
- **snake_case** for Python CLI scripts

## Database Conventions
- **Sequelize ORM** for database interactions
- **SQL queries** mixed with ORM calls
- **Manual connection pooling** in database.js
- **Validator** module used for input validation

## CLI Tool Style
- **Bash scripting** with getopts for parameter parsing
- **Python 3** for data processing scripts
- **Modular design** - separate scripts for different operations
- **Error handling** with exit codes and validation

## Documentation
- **Comprehensive README** with examples
- **Inline comments** for complex logic
- **CLAUDE.md** for AI assistant instructions
- **CLI README** separate from web app documentation

## Security Practices
- **Credentials isolation** - passwords.js gitignored
- **Input validation** in models and controllers
- **File upload restrictions** and validation
- **SQL injection prevention** via parameterized queries