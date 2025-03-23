# conf.py
import os
import sys
sys.path.insert(0, os.path.abspath('.'))

# Project information
project = 'Your Project Name'
copyright = '2025, Your Name'
author = 'Your Name'
release = '1.0.0'

# General configuration
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinxcontrib.asciidoc',  # For AsciiDoc support
]

# Document patterns
source_suffix = {
    '.rst': 'restructuredtext',
    '.adoc': 'asciidoc',
    '.asciidoc': 'asciidoc',
}

# The master document
master_doc = 'index'

# List of patterns to exclude
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Theme configuration
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# AsciiDoc settings
asciidoc_attributes = {
    'icons': 'font',
    'source-highlighter': 'rouge',
    'toc': '',
    'toclevels': '2',
}

# Setup for asciidoctor execution
def setup(app):
    # Register a custom handler for AsciiDoc files
    app.add_config_value('asciidoc_attributes', {}, 'env')