#!/usr/bin/env python
# -*- coding: utf-8 -*-

# {# pkglts, pysetup.kwds
# format setup arguments

from os import walk
from os.path import abspath, normpath
from os.path import join as pj

from setuptools import setup, find_namespace_packages


short_descr = "HydroShoot is an FSPM model to simulate gas-exchange on vine"
readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')


# find version number in src/hydroshoot/version.py
_version = {}
with open("src/openalea/hydroshoot/version.py") as fp:
    exec(fp.read(), _version)

version = _version["__version__"]

data_files = []

nb = len(normpath(abspath("src/openalea/hydroshoot_data"))) + 1


def data_rel_pth(pth):
    """ Return path relative to pkg_data
    """
    abs_pth = normpath(abspath(pth))
    return abs_pth[nb:]


for root, dnames, fnames in walk("src/openalea/hydroshoot_data"):
    for name in fnames:
        data_files.append(data_rel_pth(pj(root, name)))


setup_kwds = dict(
    name='hydroshoot',
    version=version,
    description=short_descr,
    long_description=readme + '\n\n' + history,
    author="Rami Albasha, Christian Fournier, Christophe Pradal, Eric Lebon, ",
    author_email="rami albasha at inra dot fr, @fournier-ch, @pradal, eric dot lebon at inra dot fr, ",
    url='https://github.com/openalea/hydroshoot',
    license='CeCILL-C',
    zip_safe=False,

    packages=find_namespace_packages(where='src', include=['openalea', 'openalea.*']),
    package_dir={'': 'src'},
    namespace_packages=['openalea'],

    include_package_data=True,
    package_data={'hydroshoot_data': data_files},
    install_requires=[
        ],
    tests_require=[
        "pytest",
        ],
    entry_points={},
    keywords=['FSPM', 'openalea', 'plant', 'MTG', 'Hydraulic Structure', 'gas-exchange', 'energy balance'],
    test_suite='nose.collector',
)
# #}
# change setup_kwds below before the next pkglts tag

# do not change things below
# {# pkglts, pysetup.call
setup(**setup_kwds)
# #}
