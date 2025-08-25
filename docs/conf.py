import os
import sys

project = 'JeilJungGroupCodes'
author = 'Jeil Jung Group'
extensions = [
    'sphinx.ext.mathjax',
    'sphinxfortran.fortran_domain',
    'sphinxfortran.fortran_autodoc',
]

templates_path = ['_templates']
exclude_patterns = ['_build']
html_theme = 'sphinx_rtd_theme'

# sphinx-fortran configuration
fortran_src = [
    ('jjgc', [os.path.abspath('../lanczosKuboCode_jinwoo/Src')]),
]
