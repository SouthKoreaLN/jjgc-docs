import os
import sys

project = 'JeilJungGroupCodes'
author = 'Jeil Jung Group'
extensions = [
    'breathe',
    'sphinx.ext.mathjax',
    'sphinxfortran.fortran_domain',
]

templates_path = ['_templates']
exclude_patterns = ['_build']
html_theme = 'sphinx_rtd_theme'

# Breathe configuration
breathe_projects = {
    'jjgc': os.path.abspath('_build/doxygen/xml'),
}
breathe_default_project = 'jjgc'

# Map Fortran to the Fortran domain so Breathe renders correctly
breathe_domain_by_extension = {
    'f90': 'fortran',
    'F90': 'fortran',
}
