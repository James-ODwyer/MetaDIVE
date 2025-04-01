# ==============================================================================
# MetaDIVE Sphinx Configuration File
# ==============================================================================

# ??? BUILD COMMANDS:
# ------------------------------------------------------------------------------
# ?? Build locally (for testing or file:// viewing):
# sphinx-build -b html
#
# ?? Build for GitHub Pages (relative path safe under /MetaDIVE/):
# GITHUB_PAGES=1 sphinx-build -b html docs/source docs/html
# cp -r docs/html/* docs/
# ==============================================================================

import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

# -- Project information -----------------------------------------------------
project = 'MetaDIVE'
copyright = "2025, James O'Dwyer"
author = "James O'Dwyer"
release = '1.0.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.duration',
]

templates_path = ['_templates']
exclude_patterns = []
language = 'English'

# -- Environment-specific options --------------------------------------------
on_github_pages = os.environ.get("GITHUB_PAGES", "") == "1"

if on_github_pages:
    html_baseurl = "https://james-odwyer.github.io/MetaDIVE/"
else:
    html_baseurl = "/"

# -- HTML output -------------------------------------------------------------
html_theme = 'alabaster'
html_static_path = ['_static']
