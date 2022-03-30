from setuptools import setup, find_packages
import os
import re
import sys

with open('README') as fd:
    long_description = fd.read()
    
setup(
    name='Sapphire',
    version='1.0.0',
    url='https://github.com/kcl-tscm/Sapphire',
    description='A pythonic post-processing environment for the analysis on NanoAlloys',
    author='Robert Michael Jones',
    author_email='Robert.M.Jones@kcl.ac.uk',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'ase',
        'scipy',
        'ruptures',
        'ovito',
        'sklearn'
    ],
    extras_require={'plotting': ['matplotlib', 'jupyter', 'seaborn']},
    setup_requires=['pytest-runner', 'flake8'],
    tests_require=['pytest'],
    long_description=long_description
    #entry_points={
    #    'console_scripts': ['my-command=exampleproject.example:main']
    #}
    #package_data={'exampleproject': ['data/schema.json']}
)

