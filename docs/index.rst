.. handygenome documentation master file, created by
   sphinx-quickstart on Wed Aug 16 20:47:57 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to handygenome's documentation!
=======================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Install
=======

.. code-block:: bash
    
    conda env create -n handygenome rpy2 bioconductor-dnacopy
    pip install handygenome
    
| Including all dependencies in a conda package could not be achived due to package conflicts.
| Above commands are a workaround which works well: 1) Install rpy2 and R DNAcopy within a conda environment, 2) then install other python packages with pip
    

Dependencies
============

- numpy
- pandas
- matplotlib
- biopython >= 1.80
- pysam
- pyranges
- scipy
- scikit-learn
- PyYAML
- jupyterlab
