from setuptools import setup, find_packages
import os
import re
import sys

with open('README') as fd:
    long_description = fd.read()
    
setup(
    name='Sapphire',
    version='1.0.0',
    url='https://github.com/JonesRobM/Sapphire.git',
    description='A pythonic post-processing environment for the analysis on NanoAlloys',
    author='Robert Michael Jones',
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    author_email='Robert.M.Jones@kcl.ac.uk',
    package_dir={'': 'main/Sapphire'},
    packages=find_packages(where="main/Sapphire"),
    install_requires=[
        'numpy',
        'networkx',
        'numba',
        'networkx',
        'pandas',
        'ase',
        'scipy',
        'sklearn', 
        'ruptures',
        'tensorflow',
        'pygdm2',
        'mir-flare',
        'ray',
        'scikit-learn'
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

