MESHMAGICK
==========

|Build Status| |Coverage Statu|

*MESHMAGICK* is a command line utility as well as a python module for the manipulation of meshes encountered in the hydrodynamics community.

Its primary usage is as a conversion tool between major file formats for hydrodynamic computations tools (Nemoh, Wamit, Hydrostar or Diodore) and visualization tools (stl, Tecplot, Paraview). It will be particularly useful for code to code comparisons or benchmarking.

Meshmagick also comes with several mesh manipulation capabilities: translation, rotation, scaling, clipping by a plane, symmetry, normals flipping, normals healing (making them consistent across the mesh and outgoing), cleaning (duplicate nodes merging, renumbering...).

As of the release 1.0, meshmagick provides usefull options for hydrostatics computations. It can solve for hydrostatics equilibrium for a given mass, center of gravity or both and provide the clipped mesh to be used by BEM softwares as well as the hydrostatics parameters (stiffness matrix, position of the center of buoyancy, displacement, draft...). The implementation of the algorithms for hydrostatic equilibrium resolution are iterative, highly optimized and usually converge in only a few iterations.

It also comes with an automatic meshing of lids capability, which is usefull for irregular frequencies removal by BEM softwares. This feature relies on the meshpy project which is a python wrapper around the Triangle software from Jonathan Richard Shewchuk (https://www.cs.cmu.edu/~quake/triangle.html).

Meshmagick is developped at LHEEA lab (Ecole Centrale de Nantes) in Python. It is released under the open source CeCILL licence which in practice is very similar to GPL: http://www.cecill.info/index.en.html


Basic Usage
===========

Complete this section

Installation
============

Windows
-------

You can install meshmagick by double-clicking on the install.bat file.
That's it.

UNIX
----

You can install MESHMAGICK with the following command issued:

``$ pip install .``

Updating installation with a new version
----------------------------------------

The best way to update your MESHMAGICK installation is to issue the following command in the meshmagick download directory:

``$ pip install . --upgrade``


.. |Build Status| image:: https://d-ice.githost.io/meshmagick/meshmagick/badges/release1.0/build.svg
   :target: https://d-ice.githost.io/meshmagick/meshmagick/commits/release1.0
.. |Coverage Status| image:: https://d-ice.githost.io/meshmagick/meshmagick/badges/release1.0/coverage.svg
     :target: (https://d-ice.githost.io/meshmagick/meshmagick/commits/release1.0
