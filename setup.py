"""Setup script for meshmagick."""
from __future__ import print_function
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import codecs
import os
import sys
import re

HERE = os.path.abspath(os.path.dirname(__file__))

def read(*parts):
    """Return multiple read calls to different readable objects as a single
    string."""
    # intentionally *not* adding an encoding option to open
    return codecs.open(os.path.join(HERE, *parts), 'r').read()

LONG_DESCRIPTION = read('README.rst')

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = [
            '--strict',
            '--verbose',
            '--tb=long',
            'tests']
        self.test_suite = True

    def run_tests(self):
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)

setup(
    name='meshmagick',
    version='1.0',
    url='https://', # TODO : A completer
    author='Francois Rongere -- Ecole Centrale de Nantes',
    author_email='Francois.Rongere@ec-nantes.fr',
    description="""An utility for unstructured mesh conversion and manipulation
                   for the hydrodynamic community""",
    long_description=LONG_DESCRIPTION,
    license='CeCILL-2.1',
    keywords='hydrodynamics, unstructured mesh, conversion, manipulation',
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    tests_require=['pytest', 'pytest-cov'],
    install_requires=[
        'argparse',
        'argcomplete',
        'numpy'
        ],
    cmdclass={'test': PyTest},
    entry_points={
        'console_scripts': [
            'meshmagick=meshmagick:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Development Status :: 4 - Beta',
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)',
        'Operating System :: OS Independent',
    ],
)
