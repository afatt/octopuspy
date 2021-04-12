import sys
from setuptools import setup, find_packages

__version__ = '1.0.0'

REQUIRES = [
    'coverage',
    'cycler',
    'kiwisolver',
    'matplotlib',
    'numpy',
    'pyparsing',
    'python-dateutil',
    'six'
]

setup(
    name='octopuspy',
    author='Austin Fatt',
    author_email='afatt90@gmail.com, austin.d.fatt@navy.mil',
    url='https://github.com/afatt/octopuspy',
    version=__version__,
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3.7'
    ],
    keywords='octopuspy',
    packages=find_packages(),
    install_requires=REQUIRES
)
