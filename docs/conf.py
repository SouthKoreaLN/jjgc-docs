import os
import sys

project = 'JeilJungGroupCodes'
author = 'Jeil Jung Group'
extensions = [
    'breathe',
    'sphinx.ext.mathjax',
]

templates_path = ['_templates']
exclude_patterns = ['_build']
html_theme = 'sphinx_rtd_theme'

# Breathe configuration
breathe_projects = {
    'jjgc': os.path.abspath('_build/doxygen/xml'),
}
breathe_default_project = 'jjgc'


