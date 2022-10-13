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

# import mock

# MOCK_MODULES = ['numpy', 'pandas', 'scikit-image']
# for mod_name in MOCK_MODULES:
#     sys.modules[mod_name] = mock.Mock()

# sys.path.insert(0, os.path.abspath('.'))
# If you have setup your sphinx project to use separate build and source directories, that call should instead be
sys.path.append(os.path.abspath('../..'))

# -- Project information -----------------------------------------------------

project = 'BiMaTo'
copyright = '2022, Tony Fischer (tku137)'
author = 'Tony Fischer (tku137)'

# The full version, including alpha/beta/rc tags
release = "2022.1.2"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True

# Example configuration for intersphinx: refer to the Python standard library.
# intersphinx_mapping = {'https://docs.python.org/': None}
intersphinx_mapping = {'python'    : ('http://docs.python.org/3', None),
                       'numpy'     : ('http://docs.scipy.org/doc/numpy/', None),
                       'scipy'     : ('http://docs.scipy.org/doc/scipy/reference/', None),
                       'pandas'    : ('https://pandas.pydata.org/docs/', None),
                       'matplotlib': ('http://matplotlib.sourceforge.net/', None)}

autodoc_mock_imports = [
    "numpy", "pandas", "scikit-image",
    "skimage", "skimage.filters",
    "scipy", "scipy.ndimage", "scipy.ndimage.morphology"
    ]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']