#!/usr/bin/env python

from distutils.core import setup
#from setuptools import setup

# Get version number from package
exec(open('pathfinderlib/version.py').read())

setup(
    name='PathFinder',
    version=__version__,
    description='A short wrapper for finding presence/absence/abundance of pathways from metagenomic data.',
    author='Ashwini Patil',
    author_email='patil.ashwini1091@gmail.com',
    url='https://github.com/PennChopMicrobiomeProgram',
    packages=['pathfinderlib'],
    scripts=['scripts/pathfinder.py', 'scripts/make_index_pathFinder.py'],
    install_requires=["pandas", "biopython"]
    )
