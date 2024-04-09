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
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.githubpages",
    "sphinxext.opengraph",
    "sphinx_autodoc_typehints",
    "readthedocs_ext.readthedocs",
    "sphinx_copybutton",
    "nbsphinx",
    "docutils",
]

ogp_site_url = "https://scfates.readthedocs.io/"
ogp_image = "https://scfates.readthedocs.io/en/latest/_images/scFates_logo_dark.png"

# Generate the API documentation when building
autosummary_generate = True
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [("Params", "Parameters")]

intersphinx_mapping = dict(
    python=("https://docs.python.org/3", None),
    anndata=("https://anndata.readthedocs.io/en/latest/", None),
    scanpy=("https://scanpy.readthedocs.io/en/latest/", None),
    cuml=("https://docs.rapids.ai/api/cuml/stable/", None),
    sklearn=("https://scikit-learn.org/stable/", None),
    cugraph=("https://docs.rapids.ai/api/cugraph/stable/", None),
    scipy=("https://docs.scipy.org/doc/scipy/reference/", None),
    cellrank=("https://cellrank.readthedocs.io/en/stable/", None),
    seaborn=("https://seaborn.pydata.org/", None),
)

intersphinx_disabled_domains = ["std"]

templates_path = ["_templates"]
source_suffix = [".rst", ".ipynb"]
master_doc = "index"

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

# -- Options for EPUB output
epub_show_urls = "footnote"
