Installation Instructions
=========================

Meshmagick works on every major OS (\*nix, Windows, OS/X) for both 32/64 bit systems.

The major advice to give is to rely on `Anaconda <https://www.anaconda.com/products/individual>`_ or
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ python distribution.

There is one and only one supported installation procedure that is currently and it relies on ``pip`` package installer.

After cloning the repository from GitHub and choosing the desired version (say 3.2)::

    git clone git@github.com:LHEEA/meshmagick.git
    git checkout 3.2

All you need to do from the root directory of the repository is::

    pip install -e .

Note that the ``-e`` option stands for editable and it says pip not to copy files from the repository to the python conda
install but keep the files in place. Although not mandatory, this option is very usefull as to update meshmagick to a
newer version, you only have to pull it from GiHub and it will be readily available to you, without having to reinstall.


(Optional) Run ``pytest`` if you have `pytest <http://doc.pytest.org/en/latest/>`_ installed.

