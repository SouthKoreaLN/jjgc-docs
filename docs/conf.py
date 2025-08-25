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

# sphinx-fortran configuration (domain only, no autodoc)
fortran_src = [
    ('jjgc', [os.path.abspath('../lanczosKuboCode_jinwoo/Src')]),
]
fortran_ext = ['f90', 'F90']
