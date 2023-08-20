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

You can install handygenome by running the commands below,

.. code-block:: bash
    
    conda env create -n handygenome rpy2 bioconductor-dnacopy
    pip install handygenome

| which 1) makes a conda environment to install rpy2 using R with DNAcopy package, 2) then installs remaining python dependencies with pip.
|   

    Including all dependencies in a single conda package could not be achived due to package conflicts.
    

Python Dependencies
===================

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
