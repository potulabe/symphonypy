# Configuration file for the Sphinx documentation builder.

# -- Project information

project = "Symphonypy"
copyright = "2023, Petrova, Isaev"
authors = "Petrova, Isaev"

release = "0.2.1"
version = "0.2.1"

# -- General configuration ------------------------------------------------

needs_sphinx = "2.0"

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "nbsphinx",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}
intersphinx_disabled_domains = ["std"]

templates_path = ["_templates"]

# -- Options for HTML output

html_theme = "sphinx_rtd_theme"

# -- Options for EPUB output
epub_show_urls = "footnote"

# -- Retrieve notebooks ------------------------------------------------

from urllib.request import urlretrieve

notebooks_url = "https://github.com/potulabe/symphonypy/raw/main/notebooks/"
notebooks = [
    "Symphonypy_precomputed.ipynb",
    "Symphonypy_simple_tutorial.ipynb",
    "Symphonypy_without_harmony_tutorial.ipynb",
]
for nb in notebooks:
    try:
        urlretrieve(notebooks_url + nb, nb)
    except:
        pass
