import sys
from setuptools import setup, find_packages

__version__ = '1.0.0'

REQUIRES = [
    'cycler==0.10.0',
    'kiwisolver==1.3.1',
    'matplotlib==3.1.0',
    'numpy==1.20.1',
    'pyparsing==2.4.7',
    'python-dateutil==2.8.1',
    'six==1.15.0'
]

setup(
    name='octopuspy',
    author='Austin Fatt',
    author_email='austin.d.fatt@navy.mil',
    url='',
    version=__version__,
    license='',
    classifiers=[
        'Programming Language :: Python :: 3.7'
    ],
    keywords='octopuspy',
    packages=find_packages(),
    install_requires=REQUIRES
)
