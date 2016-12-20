Meshmagick
==========

**Meshmagick** is a command line utility as well as a python module for the manipulation of meshes encountered in the
hydrodynamics community.

Its primary goal was to be a conversion tool between major file formats for hydrodynamic computations tools (Nemoh,
Wamit, Hydrostar or Diodore) and visualization tools (stl, Tecplot, Paraview). It will be particularly useful for code
to code comparisons or benchmarking.

**Meshmagick** also comes with several mesh manipulation capabilities: translation, rotation, scaling, clipping by a
plane, symmetry, normals flipping, normals healing (making them consistent across the mesh and outgoing), cleaning
(duplicate nodes merging...).

As of the release 1.0, **meshmagick** provides useful options for hydrostatics computations. It can solve for
hydrostatics equilibrium for a given mass, center of gravity or both and provide the clipped mesh to be used by BEM
software as well as the hydrostatics parameters (stiffness matrix, position of the center of buoyancy, displacement,
draft...). Inertial properties of meshes may also be computed, based on assumptions.

**Meshmagick** is primarily a command line utility for everyday hydrodynamicists. However, it also comes with a
package that can be imported in a python script and give the full access to the command line options, programmatically.

.. note::
    **Meshmagick is the property of Ecole Centrale de Nantes and is maintained by François Rongère <francois
    .rongere@ec-nantes.fr>**. It is released under the GNU GPL v3 open source licence (see LICENCE file).

Getting Meshmagick
------------------

.. warning::
    Meshmagick is written in Python 2.7 only. So please ensure that you have a Python 2.7 distribution.

Meshmagick is packaged for every major OS (*nix, Windows, OS/X) for both 32/64 bit systems.

Installing with conda
~~~~~~~~~~~~~~~~~~~~~

This is the prefered method to get your own packaged copy of Meshmagick as `conda <http://conda.pydata.org/docs/>`_
is known to deal with every dependencies of Meshmagick. For this, you will need to install
`Anaconda <https://www.continuum.io/downloads>`_ Python distribution or simply
`Miniconda <http://conda.pydata.org/miniconda.html>`_ for your architecture which are shipped with the conda package
manager. If you are not an everyday Python developer, maybe you will prefer installing Miniconda which is far mush
lighter than Anaconda.

Once conda installed, just try this::

    conda install -c frongere meshmagick

It will download the Meshmagick package corresponding to your architecture. If everything goes well, it will install
for you the dependencies. You may now test the install by::

    meshmagick -h

which should now show you the embedded help of the command line tool.

Updating to the latest meshmagick version is obtained by::

    conda update meshmagick

Installing with pip
~~~~~~~~~~~~~~~~~~~

For those wanting to use pip, just do::

    pip install meshmagick

It will install Meshmagick from Pypi.

Update to the newest version is achieved by::

    pip install meshmagick --upgrade

.. note::
    This is not the preferred method as Meshmagick depends on vtk for its visualization features which is sadly not
    available on Pypi. You will then have to deal with this dependency by yourself.

Installing from source
~~~~~~~~~~~~~~~~~~~~~~

This can be done by checking out the source files from the Git source code repository. This option is mainly for
those wanting to develop into Meshmagick.

1. Clone the meshmagick repository::

    git clone https://github.com/LHEEA/meshmagick.git

2. Change directory to meshmagick

3. Run::

    python setup.py install

4. (Optional) Run ``pytest`` if you have `pytest <http://doc.pytest.org/en/latest/>`_ installed.

