# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import sys
from pathlib import Path
sys.path.insert(0, str(Path('../..').resolve()))

project = 'GridCellsCond'
copyright = '2024, Inayath Shaikh'
author = 'Inayath Shaikh'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autosummary',
                'sphinx.ext.duration',
                'sphinx.ext.napoleon',
                "sphinx_inline_tabs",
                'sphinx.ext.imgconverter'
    ]

templates_path = ['_templates']
exclude_patterns = ['s_sim_setup.py']
napoleon_use_param = False
html_title = "GridCellsCond"
html_theme = "furo"
html_static_path = ['_static']
