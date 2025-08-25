# Complete Working Documentation Setup for Fortran 90

This document explains the complete working configuration for generating documentation from Fortran 90 files using Doxygen + Breathe + Sphinx + ReadTheDocs.

## Overview

The setup successfully renders Fortran 90 functions with full Doxygen documentation in ReadTheDocs using:
- **Doxygen**: Parses Fortran 90 files and generates XML
- **Breathe**: Bridges Doxygen XML to Sphinx
- **Sphinx**: Renders the final HTML documentation
- **ReadTheDocs**: Hosts and builds the documentation

## Why This Setup Works

### The Problem
Doxygen has fundamental limitations with Fortran 90 parsing. Even with proper `!>` Doxygen headers, it often fails to extract individual function documentation, resulting in blank pages.

### The Solution
1. **Restrict Doxygen input** to only the specific file(s) you want documented
2. **Force extraction** with `EXTRACT_ALL=YES`
3. **Use proper namespace rendering** instead of file indexing
4. **Bridge through Breathe** for reliable Sphinx integration

## Complete Working Configuration

### 1. Doxyfile (Key Settings)

```bash
PROJECT_NAME           = "JeilJungGroupCodes"
PROJECT_NUMBER         = 1.0
OUTPUT_DIRECTORY       = docs/_build/doxygen
GENERATE_XML           = YES
XML_OUTPUT             = xml
GENERATE_HTML          = NO
RECURSIVE              = YES

# CRITICAL: Restrict input to avoid noise
INPUT                  = lanczosKuboCode_jinwoo/Src/ham.F90

FILE_PATTERNS          = *.f90 *.F90
EXTENSION_MAPPING      = F90=Fortran f90=Fortran
OPTIMIZE_OUTPUT_FOR_C  = NO
OPTIMIZE_FOR_FORTRAN   = YES

# CRITICAL: Force extraction of all documented items
EXTRACT_ALL            = YES
EXTRACT_PRIVATE        = NO
QUIET                  = NO
WARN_IF_UNDOCUMENTED   = YES
EXTRACT_LOCAL_CLASSES  = YES
EXTRACT_LOCAL_METHODS  = YES
EXTRACT_ANON_NSPACES   = YES
```

**Why These Settings Matter:**
- `INPUT = lanczosKuboCode_jinwoo/Src/ham.F90`: Focuses Doxygen on one file, avoiding noise from other files
- `EXTRACT_ALL = YES`: Forces Doxygen to extract all documented items, even if it's unsure about them
- `OPTIMIZE_FOR_FORTRAN = YES`: Enables Fortran-specific parsing optimizations

### 2. Fortran 90 Doxygen Syntax (Working Format)

```fortran
!> @brief Brief description of what the function does.
!! @param[in]  param1   Description of input parameter
!! @param[out] param2   Description of output parameter
!! @param[in]  param3   Description of another parameter
!! @return Description of return value (if applicable)
!! @note Additional notes about the function
!! @see relatedFunction
subroutine functionName(param1, param2, param3)
    ! Function implementation
end subroutine functionName
```

**Key Points:**
- Use `!>` for the main description
- Use `!!` for additional lines
- Always document parameters with `@param[in/out]`
- Include `@brief` for short descriptions

### 3. Sphinx Configuration (docs/conf.py)

```python
import os
import sys

project = 'JeilJungGroupCodes'
author = 'Jeil Jung Group'
extensions = [
    'breathe',           # CRITICAL: Bridges Doxygen to Sphinx
    'sphinx.ext.mathjax',
]

templates_path = ['_templates']
exclude_patterns = ['_build']
html_theme = 'sphinx_rtd_theme'

# Breathe configuration - CRITICAL for Doxygen integration
breathe_projects = {
    'jjgc': os.path.abspath('_build/doxygen/xml'),  # Points to Doxygen XML output
}
breathe_default_project = 'jjgc'
```

### 4. API Documentation (docs/api.rst)

```rst
API Reference
=============

ham module
----------

.. doxygennamespace:: ham
   :project: jjgc
   :members:
   :undoc-members:
   :no-link:
```

**Why This Works:**
- `.. doxygennamespace:: ham` renders the entire namespace with all its members
- `:members:` shows all documented members
- `:undoc-members:` includes even undocumented members (useful for debugging)

### 5. Index Page (docs/index.rst)

```rst
JeilJungGroupCodes documentation
================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

API Reference
-------------

.. toctree::
   :maxdepth: 2

   api
```

**Critical Fix:**
- Changed from `.. doxygenindex::` (shows file lists) to including `api.rst`
- This ensures function documentation is rendered instead of blank pages

### 6. ReadTheDocs Configuration (.readthedocs.yaml)

```yaml
version: 2
build:
  os: ubuntu-22.04
  tools:
    python: "3.11"
  apt_packages:
    - doxygen  # CRITICAL: Install Doxygen
  jobs:
    pre_build:
      - doxygen Doxyfile  # CRITICAL: Run Doxygen before Sphinx

sphinx:
  configuration: docs/conf.py

python:
  install:
    - requirements: docs/requirements.txt

formats:
  - htmlzip
```

### 7. Python Requirements (docs/requirements.txt)

```
Sphinx==7.3.7
breathe==4.35.0
sphinx-rtd-theme==2.0.0
sphinx-fortran==1.1.1
six==1.16.0
numpy==2.1.2
```

## Troubleshooting Guide

### Problem: Blank Pages in Documentation
**Symptoms:** Functions are documented in Fortran files but don't appear in rendered docs

**Solutions (in order):**
1. **Check Doxygen XML output**: Look for `docs/_build/doxygen/xml/` directory
2. **Verify function extraction**: Search XML for function names (they'll be lowercase)
3. **Check INPUT restriction**: Ensure Doxygen only processes the files you want
4. **Verify EXTRACT_ALL=YES**: This forces extraction of all documented items
5. **Check namespace rendering**: Use `.. doxygennamespace::` not `.. doxygenindex::`

### Problem: Doxygen Not Finding Functions
**Symptoms:** XML output is empty or missing function definitions

**Solutions:**
1. **Restrict INPUT**: Use `INPUT = path/to/specific/file.F90` instead of entire directories
2. **Enable EXTRACT_ALL**: Set `EXTRACT_ALL = YES` in Doxyfile
3. **Check file patterns**: Ensure `FILE_PATTERNS = *.f90 *.F90` includes your files
4. **Verify syntax**: Check that Doxygen headers use `!>` and `!!` correctly

### Problem: Breathe Can't Read XML
**Symptoms:** Sphinx build fails with Breathe errors

**Solutions:**
1. **Check XML path**: Ensure `breathe_projects` points to correct XML directory
2. **Verify Doxygen ran**: Check that `docs/_build/doxygen/xml/` exists and has content
3. **Check XML format**: Ensure Doxygen generated valid XML (not just HTML)

## Common Pitfalls

### ❌ Don't Do This
- Use `INPUT = lanczosKuboCode_jinwoo/Src` (too broad, creates noise)
- Use `.. doxygenindex::` (shows files, not functions)
- Set `QUIET = YES` (hides important warnings)
- Skip `EXTRACT_ALL = YES` (may miss documented items)

### ✅ Do This Instead
- Use `INPUT = lanczosKuboCode_jinwoo/Src/ham.F90` (focused, clean)
- Use `.. doxygennamespace:: ham` (shows functions with docs)
- Set `QUIET = NO` (shows warnings for debugging)
- Always use `EXTRACT_ALL = YES` (ensures extraction)

## Testing the Setup

### Local Testing
```bash
# 1. Generate Doxygen XML
doxygen Doxyfile

# 2. Check XML output
ls docs/_build/doxygen/xml/

# 3. Search for functions in XML
grep -r "functionname" docs/_build/doxygen/xml/

# 4. Build Sphinx locally (if dependencies installed)
cd docs && python -m sphinx -b html . _build/html
```

### ReadTheDocs Testing
1. Push changes to your docs repository
2. Check ReadTheDocs build logs for errors
3. Verify that functions appear in the rendered documentation
4. Check that parameter documentation is visible

## File Structure

```
your-docs-repo/
├── .readthedocs.yaml          # ReadTheDocs configuration
├── Doxyfile                   # Doxygen configuration
├── docs/
│   ├── conf.py               # Sphinx configuration
│   ├── index.rst             # Main index page
│   ├── api.rst               # API documentation
│   ├── requirements.txt      # Python dependencies
│   └── DOCUMENTATION_SETUP.md # This file
└── lanczosKuboCode_jinwoo/
    └── Src/
        └── ham.F90           # Your documented Fortran file
```

## Summary

The key to success is:
1. **Focus Doxygen** on specific files only
2. **Force extraction** with `EXTRACT_ALL=YES`
3. **Use proper namespace rendering** in Sphinx
4. **Bridge through Breathe** for reliable integration

This setup successfully renders Fortran 90 functions with full Doxygen documentation in ReadTheDocs, overcoming Doxygen's inherent limitations with Fortran parsing.

## Version History

- **Initial Setup**: Basic Doxygen + Breathe + Sphinx
- **First Fix**: Restricted INPUT to specific files
- **Second Fix**: Enabled EXTRACT_ALL=YES
- **Final Fix**: Changed from doxygenindex to doxygennamespace rendering
- **Current**: Working configuration with comprehensive documentation

---

*This document should be updated whenever the configuration changes or new troubleshooting steps are discovered.*
