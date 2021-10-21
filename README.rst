========================
HydroShoot
========================

.. {


.. image:: https://travis-ci.org/openalea/hydroshoot.svg?branch=master
    :alt: Travis build status
    :target: https://travis-ci.org/openalea/hydroshoot

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
    ``conda create -n MyEnv -c fredboudon -c conda-forge openalea.plantgl alinea.caribu alinea.astk python=3.8``

    ``conda activate MyEnv``

    ``conda install -c pvlib pvlib``

    ``conda install -c conda-forge sphinx sphinx-gallery sphinx_rtd_theme``

    ``conda install pandas numpy scipy sympy jsonschema matplotlib pytz``

    ``cd ~/DIRECTORY-WHERE-HYDROSHOOT-WILL-BE-INSTALLED/``

    ``git clone https://github.com/openalea/hydroshoot.git``

    (Up to you to use SSH or GitHub CLI inseated)

    ``cd hydroshoot/``

    ``python setup.py develop``

You're done !


