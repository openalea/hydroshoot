========================
hydroshoot
========================

.. {# pkglts, doc


.. image:: https://travis-ci.org/Rami Albasha/hydroshoot.svg?branch=master
    :alt: Travis build status
    :target: https://travis-ci.org/Rami Albasha/hydroshoot

.. #}

Installation
------------

    conda install hydroshoot -c default -c openalea -c conda-forge


**Tempo**
HydroShoot is a functional-structural plant modelling package. It is built using the MTG (Multiscale Tree Graph) central data structure (Pradal et al., 2008).

Hydroshoot is composed of 4 generic modules used to simulate:
	- irradiance interception (*irradiance*)
	- xylem water transport (*hydraulic*)
	- leaves energy budget (*energy*)
	- leaves gas exchange fluxes (*exchange*).

Hydroshoot provides a grapevine-specific module (*architecture*) which builds plant shoot structure for potted of trained grapevines.


