# Arotake

**Arotake** short description.

## Public API

## Overview

## Usage

## Development

### Code Quality and Testing Commands

- `make fmt` - Runs ruff format, which automatically reformats Python files according to the style rules in `pyproject.toml`
- `make lint` - Runs ruff check --fix, which lints the code (checks for style errors, bugs, outdated patterns, etc.) and auto-fixes what it can.
- `make type` - Runs mypy, the static type checker, using the strictness settings from `pyproject.toml`. Mypy is a static type checker for Python, a dynamically typed language. Because static analysis cannot account for all dynamic runtime behaviors, mypy may report false positives which do no reflect actual runtime issues. The usefulness of mypy is therefore limited, unless the developer compensates with extra work for the choices that were made when Python was originally designed.
- `make test` - Runs pytest with coverage reporting (configured in `pyproject.toml`).

## Author

Jan Kazil - jan.kazil.dev@gmail.com - [jankazil.com](https://jankazil.com)

## License

BSD-3-Clause

