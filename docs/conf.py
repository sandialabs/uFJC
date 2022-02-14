# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath("../../"))

# -- Project information -----------------------------------------------------

project = 'uFJC'
copyright = '2022, Sandia National Laboratories'
author = 'Michael R. Buche, Scott J. Grutzik'

# The full version, including alpha/beta/rc tags
release = '1.0.0'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'matplotlib.sphinxext.plot_directive',
    'sphinxcontrib.bibtex'
]

# bibtex files for citing in docs and bib style
bibtex_bibfiles = ['main.bib']
bibtex_default_style = 'plain'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# Show __init__ methods in the documentation
autodoc_default_options = {
    'special-members': '__init__',
}

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'navigation_depth': 8
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# For LaTeX math
latex_engine = 'xelatex'

# Opt not to show HTML links when including matplotlib plots
plot_html_show_formats = False
plot_html_show_source_link = False

# But do show code blocks used to make the matplotlib plots
plot_include_source = True

# Do not add module name(s) before classes, methods
add_module_names = False
