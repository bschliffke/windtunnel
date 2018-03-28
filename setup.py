# -*- coding: utf-8 -*-
import sys
from distutils.core import setup
from os.path import dirname, join
from setuptools import find_packages


if not sys.version_info >= (3, 5, 1):
    sys.exit('Only supports Python version >=3.5.1.\n'
             'Current version is {}'.format(sys.version))

__version__ = open(join(dirname(__file__), 'windtunnel', 'VERSION')).read().strip()

setup(
    name='windtunnel',
    author='Benyamin Schliffke',
    author_email='benny.schliffke@gmail.com',
    url='https://github.com/bschliffke/windtunnel',
    download_url='https://github.com/bschliffke/windtunnel/' + __version__,
    version=__version__,
    packages=find_packages(),
    license='MIT',
    description='Python package for use with BSA software output.',
    classifiers=[
        # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    include_package_data=True,
    install_requires=[
        'fnmatch',
        'logging',
        'math',
        'matplotlib>=1.5.1',
        'numpy>=1.10.4',
        'os',
        'pandas',
        'scipy',
        'skimage',
    ],
)