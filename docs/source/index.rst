Symphonypy — port of Symphony algorithm of single-cell reference atlas mapping to Python
========================================================================================

Symphonypy is a Python port of the original R Symphony package for label transfer.

In the case of usage in your work, please cite the original Symphony paper::

    Kang JB, Nathan A, Weinand K, Zhang F, Millard N, Rumker L, Moody DB, Korsunsky I, Raychaudhuri S
    Efficient and precise single-cell reference atlas mapping with Symphony
    Nat Commun. 2021 Oct 7;12(1):5890. doi: 10.1038/s41467-021-25957-x

The package inculdes slightly modified version of ``sce.pp.harmony_integrate()`` and ``sc.tl.ingest()`` for 
better experience and gives an interface for t-SNE embedding coordinates transfer via ``openTSNE`` Python package.

.. toctree::
   :caption: Main
   :maxdepth: 2
   :hidden:

   usage
   api

.. toctree::
   :caption: Walkthroughs
   :maxdepth: 2
   :hidden:

   Symphonypy_precomputed
   Symphonypy_simple_tutorial
   Symphonypy_without_harmony_tutorial
