Installation Instructions
=========================

Meshmagick works on every major OS (\*nix, Windows, OS/X) for both 32/64 bit systems.

.. warning::

    From version 2.0, Meshmagick is now written in Python 3. Compatibility with Python 2.7 is no more ensured. If special
    needs for this compatibility, please contact the developers. Thanks.


Installing with conda
---------------------

This is the prefered method to get your own packaged copy of Meshmagick as `conda <http://conda.pydata.org/docs/>`_
is known to deal with every dependencies of Meshmagick. For this, you will need to install
`Anaconda <https://www.continuum.io/downloads>`_ Python distribution or simply
`Miniconda <http://conda.pydata.org/miniconda.html>`_ for your architecture which are shipped with the conda package
manager. If you are not an everyday Python developer, maybe you will prefer installing Miniconda which is far mush
lighter than Anaconda.

.. warning::

    **For Windows users**, current version of conda (4.3.*) appears to be buggy with respect to conda install of
    meshmagick. One workaround is to downgrade conda to version 4.2.* by typing::

        conda install conda=4.2

Once conda installed, just try this::

    conda install -c frongere meshmagick

It will download the Meshmagick package corresponding to your architecture. If everything goes well, it will install
for you the dependencies. You may now test the install by::

    meshmagick -h

which should now show you the embedded help of the command line tool.

Updating to the latest meshmagick version is obtained by::

    conda update meshmagick

.. warning::

    This method for installation is quick. However, we do not systematically package Meshmagick for conda. We apologize
    for that. The advice is to install Meshmagick from sources with pip using the edit mode -e (see below).


Installing with pip
-------------------

For those wanting to use pip, just do::

    pip install meshmagick

It will install Meshmagick from Pypi.

Update to the newest version is achieved by::

    pip install meshmagick --upgrade

.. note::
    This is not the preferred method as Meshmagick depends on vtk for its visualization features which is sadly not
    available on Pypi. You will then have to deal with this dependency by yourself.

    * **On Linux** : You may use the current package manager of your system to install vtk
    * **On Mac** : You may follow the instructions given in *bemio* documentation:
      https://wec-sim.github.io/bemio/installing.html#installing-vtk-with-python-bindings
    * **On Windows** : We had success in installing unofficial pakage from here:
      http://www.lfd.uci.edu/~gohlke/pythonlibs/

Installing from source
----------------------

This can be done by checking out the source files from the Git source code repository. This option is mainly for
those wanting to develop into Meshmagick.

1. Clone the meshmagick repository::

    git clone https://github.com/LHEEA/meshmagick.git

2. Change directory to meshmagick

3. Install:

    * If you want to install from source for just using meshmagick, it can be achieved by::

        python setup.py install

    * If you want to install from soruce for development purposes, you should better install meshmagick in
      development mode so that your code modifications are directly taken into account with respect to your install.
      It is achieved by the command::

        pip install -e .

4. (Optional) Run ``pytest`` if you have `pytest <http://doc.pytest.org/en/latest/>`_ installed.

