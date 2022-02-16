import os
import sys
sys.path.insert(0, os.path.abspath("../../"))

release = '1.0.0'
project = 'uFJC'
author = 'Michael R. Buche, Scott J. Grutzik'
copyright = '2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS). \
    Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.'

templates_path = ['_templates']
html_static_path = ['_static']
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'navigation_depth': 8
}
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'matplotlib.sphinxext.plot_directive',
    'sphinxcontrib.bibtex'
]
latex_engine = 'xelatex'
bibtex_bibfiles = ['main.bib']
bibtex_default_style = 'plain'
autodoc_default_options = {
    'special-members': '__init__',
}
plot_html_show_formats = False
plot_html_show_source_link = False
plot_include_source = True
add_module_names = False
