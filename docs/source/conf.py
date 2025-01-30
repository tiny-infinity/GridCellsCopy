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
                # 'sphinx.ext.autodoc',
                # 'sphinx_rtd_theme',
                'sphinx.ext.napoleon',
    ]

templates_path = ['_templates']
exclude_patterns = ['s_sim_setup.py']
# napoleon_google_docstring  =False
napoleon_use_param = False
# napoleon_use_rtype = False

# autodoc_mock_imports = [
#     's_sim_setup'
# ]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = "sphinx_rtd_theme"
html_theme = "furo"
# html_theme="furo"
html_static_path = ['_static']
