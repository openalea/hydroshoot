========================
HydroShoot
========================

.. {


.. image:: https://github.com/openalea/hydroshoot/actions/workflows/conda-package-build.yml/badge.svg
    :alt: CI status
    :target: https://github.com/openalea/hydroshoot/actions/workflows/conda-package-build.yml
    
.. image:: https://readthedocs.org/projects/hydroshoot/badge/?version=latest
    :target: https://hydroshoot.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
    
.. image:: https://anaconda.org/openalea3/openalea.hydroshoot/badges/version.svg   
    :target: https://anaconda.org/openalea3/openalea.hydroshoot

.. }


Brief Description
-----------------

HydroShoot is a functional-structural plant modelling package. 

Hydroshoot is composed of 3 generic modules used to simulate:
    - xylem water transport (*hydraulic*)
    - leaves energy budget (*energy*)
    - leaves gas exchange fluxes (*exchange*).

Hydroshoot provides a grapevine-specific module (*architecture*) which builds plant shoot structure for potted of trained grapevines.



Installation
------------

HydroShoot must temporarily be installed manually by cloning this repo locally then installing.
An installation procedure using conda will be possible in the future.

You can follow the following steps for installation:
     ``mamba create -n hydroshoot -c openalea3 -c conda-forge openalea.hydroshoot``
 
     ``mamba activate hydroshoot``


You're done !


Documentation
-------------

https://hydroshoot.readthedocs.io/

Citation
--------

R Albasha, C Fournier, C Pradal, M Chelle, J A Prieto, G Louarn, T Simonneau, E Lebon, HydroShoot: a functional-structural plant model for simulating hydraulic structure, gas and energy exchange dynamics of complex plant canopies under water deficit—application to grapevine (Vitis vinifera), in silico Plants, Volume 1, Issue 1, 2019, diz007, https://doi.org/10.1093/insilicoplants/diz007
