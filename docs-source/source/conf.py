# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'CRISPY'
copyright = '2025, L. Brémaud, J. Girardot'
author = 'L. Brémaud, J. Girardot'
release = '0.2.0'
version = '0.2'

html_logo = "_static/logo.png"
#html_favicon = "_static/favicon.ico"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
   'sphinx.ext.duration',
   'sphinx.ext.doctest',
   'sphinx.ext.autodoc',
   'sphinx.ext.autosummary',
   'sphinx.ext.napoleon',
   'sphinx.ext.viewcode',
   "sphinx.ext.intersphinx",
   'sphinx.ext.githubpages',
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
# Available themes: furo, sphinx_rtd_theme, sphinx_book_theme, pydata_sphinx_theme, alabaster

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = [
    'custom.css',
]

autoclass_content = 'both'
autodoc_member_order = 'bysource'

html_theme_options = {
    'collapse_navigation': True,
    'navigation_depth': 3,
    'titles_only': False,
    'flyout_display': 'attached',
    'version_selector': True
}

pygments_style = 'sphinx'  # Options: 'default', 'sphinx', 'friendly', 'monokai', etc.

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pyvista': ('https://pyvista.org/', None),
}


