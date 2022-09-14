Meshmagick
==========

**NEW IN v3.3**: It is now possible to load last .msh file format version (4.1) of GMSH generated mesh files.
Compatibility with format version 2.1 is still ensured.


**Meshmagick** is a command line utility as well as a python module for the manipulation of meshes encountered in the
hydrodynamics community.

Its primary goal was to be a conversion tool between major file formats for hydrodynamic computations tools (`Nemoh
<https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp>`_, `Wamit <http://www.wamit.com/>`_, `Hydrostar
<http://www.veristar.com/portal/veristarinfo/detail/software/Seakeeping%20and%20Mooring%20Analysis/HYDROSTAR/Hydros>`_
or `Diodore <http://www.principia.fr/expertise-fields-software-products-diodore-132.html>`_) and visualization tools
(stl, `Tecplot <http://www.tecplot.com/>`_, `Paraview <http://www.paraview.org/>`_). It will be particularly useful for
code to code comparisons or benchmarking.

**Meshmagick** also comes with several mesh manipulation capabilities: translation, rotation, scaling, clipping by a
plane, symmetry, normals flipping, normals healing (making them consistent across the mesh and outgoing), cleaning
(duplicate nodes merging...).

**Meshmagick** provides useful options for hydrostatics computations. It can solve for
hydrostatics equilibrium for a given mass, center of gravity or both and provide the clipped mesh to be used by BEM
software as well as the hydrostatics parameters (stiffness matrix, position of the center of buoyancy, displacement,
draft...). Inertial properties of meshes may also be computed, based on assumptions.

**Meshmagick** is primarily a command line utility for everyday hydrodynamicists. However, it also comes with a
package that can be imported in a python script and give the full access to the command line options, programatically.

.. note::
    **Meshmagick** is the property of **Ecole Centrale de Nantes** and is maintained by François Rongère <francois
    .rongere@dice-engineering>. It is released under the **GNU GPLv3** open source licence (see LICENCE file).

GitHub Repository
-----------------

https://github.com/LHEEA/meshmagick

Documentation
-------------

https://lheea.github.io/meshmagick


Getting Meshmagick
------------------

Getting the latest version::

    pip install https://github.com/LHEEA/meshmagick/archive/master.zip
