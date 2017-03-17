"""Setup script for meshmagick."""
from __future__ import print_function
from setuptools import setup, find_packages
import codecs
import os

HERE = os.path.abspath(os.path.dirname(__file__))


def read(*parts):
    """Return multiple read calls to different readable objects as a single
    string."""
    # intentionally *not* adding an encoding option to open
    return codecs.open(os.path.join(HERE, *parts), 'r').read()

LONG_DESCRIPTION = read('README.rst')

setup(
    name='meshmagick',
    version='1.0.6',
    url='https://github.com/LHEEA/meshmagick',
    author='Francois Rongere -- Ecole Centrale de Nantes',
    author_email='Francois.Rongere@ec-nantes.fr',
    description="""A command line tool to manipulate hydrodynamics meshes""",
    long_description=LONG_DESCRIPTION,
    license='GPLv3',
    keywords='hydrodynamics, unstructured mesh, conversion, manipulation',
    packages=find_packages(exclude=['contrib', 'doc', 'tests*']),
    # setup_requires=['pytest-runner'],
    # tests_require=['pytest', 'pytest-cov'],
    install_requires=['numpy', 'argcomplete'],
    entry_points={
        'console_scripts': [
            'meshmagick=meshmagick:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Development Status :: 5 - Production/Stable',
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
)
