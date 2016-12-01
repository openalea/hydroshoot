#!/usr/bin/env python
# -*- coding: utf-8 -*-

# {# pkglts, pysetup.kwds
# format setup arguments

from os import walk
from os.path import abspath, normpath
from os.path import join as pj

from setuptools import setup, find_packages


short_descr = "HydroShoot is an FSPM model to simulate gas-exchange on vine"
readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')


# find version number in src/hydroshoot/version.py
version = {}
with open("src/hydroshoot/version.py") as fp:
    exec(fp.read(), version)


data_files = []

nb = len(normpath(abspath("src/hydroshoot_data"))) + 1


def data_rel_pth(pth):
    """ Return path relative to pkg_data
    """
    abs_pth = normpath(abspath(pth))
    return abs_pth[nb:]


for root, dnames, fnames in walk("src/hydroshoot_data"):
    for name in fnames:
        data_files.append(data_rel_pth(pj(root, name)))


setup_kwds = dict(
    name='hydroshoot',
    version=version["__version__"],
    description=short_descr,
    long_description=readme + '\n\n' + history,
    author="Rami Albasha, Christian Fournier, Christophe Pradal, Eric Lebon, ",
    author_email="rami albasha at inra dot fr, @fournier-ch, @pradal, eric dot lebon at inra dot fr, ",
    url='https://github.com/openalea-incubator/hydroshoot',
    license='cecill-c',
    zip_safe=False,

    packages=find_packages('src'),
    package_dir={'': 'src'},
    
    include_package_data=True,
    package_data={'hydroshoot_data': data_files},
    install_requires=[
        ],
    tests_require=[
        "mock",
        "nose",
        ],
    entry_points={},
    keywords='',
    test_suite='nose.collector',
)
# #}
# change setup_kwds below before the next pkglts tag

# do not change things below
# {# pkglts, pysetup.call
setup(**setup_kwds)
# #}
