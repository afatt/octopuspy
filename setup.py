import sys
from setuptools import setup, find_packages

__version__ = '1.0.2'

readme = 'README.md'
long_description = open(readme).read()

REQUIRES = [
    'coverage==4.3.4',
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
    description='Python package for manipulating Octopus output files',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/afatt/octopuspy',
    download_url='https://github.com/afatt/octopuspy/archive/refs/tags/{}.tar.gz'.format(__version__),
    version=__version__,
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3.7'
    ],
    keywords='octopuspy',
    packages=find_packages(),
    install_requires=REQUIRES,
    python_requires='>=3.7'
)
