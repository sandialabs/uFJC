####
uFJC
####

|build| |docs| |codecov| |coveralls| |pyversions| |pypi| |conda| |docker| |license| |zenodo|

The Python package for the uFJC single-chain model.

************
Installation
************

The package can be installed using ``pip`` via the `Python Package Index <https://pypi.org/project/ufjc/>`_ (PyPI),

::

    pip install ufjc

or using ``conda`` via the ``mrbuche`` channel on `Anaconda <https://anaconda.org/mrbuche/ufjc>`_,

::

    conda install --channel mrbuche ufjc
    
Alternatively, a branch can be directly installed using

::

    pip install git+https://github.com/sandialabs/ufjc.git@<branch-name>

or after cloning a branch and executing ``python setup.py install``.
There are also `Docker images <https://hub.docker.com/r/mrbuche/ufjc>`_ available for use.
In all of these cases, a valid installation can be tested by running

::

    python -m ufjc.tests

***********
Information
***********

- `Contributing <https://ufjc.readthedocs.io/en/latest/CONTRIBUTING.html>`__
- `Documentation <https://ufjc.readthedocs.io/en/latest/>`__
- `Examples <https://ufjc.readthedocs.io/en/latest/ufjc.examples.html>`__
- `License <https://github.com/sandialabs/ufjc/blob/main/LICENSE>`__
- `Releases <https://github.com/sandialabs/ufjc/releases>`__
- `Repository <https://github.com/sandialabs/ufjc>`__
- `Tutorial <https://ufjc.readthedocs.io/en/latest/TUTORIAL.html>`__

********
Citation
********

\M. R. Buche and S. J. Grutzik, ``uFJC``: the Python package for the uFJC single-chain model, `Zenodo (2022) <https://doi.org/10.5281/zenodo.6114263>`_.

\M. R. Buche, M. N. Silberstein, and S. J. Grutzik, Freely jointed chain models with extensible links, `Physical Review E 106, 024502 (2022) <https://doi.org/10.1103/PhysRevE.106.024502>`_.

*********
Copyright
*********

Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.

..
    Badges ========================================================================

.. |docs| image:: https://img.shields.io/readthedocs/ufjc?logo=readthedocs&label=Read%20the%20Docs
    :target: https://ufjc.readthedocs.io/en/latest/

.. |build| image:: https://img.shields.io/github/workflow/status/sandialabs/ufjc/main?label=GitHub&logo=github
    :target: https://github.com/sandialabs/ufjc

.. |coveralls| image:: https://img.shields.io/coveralls/github/sandialabs/ufjc?logo=coveralls&label=Coveralls
    :target: https://coveralls.io/github/sandialabs/ufjc?branch=main

.. |codecov| image:: https://img.shields.io/codecov/c/github/sandialabs/ufjc?label=Codecov&logo=codecov
    :target: https://codecov.io/gh/sandialabs/ufjc

.. |pyversions| image:: https://img.shields.io/pypi/pyversions/ufjc.svg?logo=python&logoColor=FBE072&color=4B8BBE&label=Python
    :target: https://pypi.org/project/ufjc/

.. |pypi| image:: https://img.shields.io/pypi/v/ufjc?logo=pypi&logoColor=FBE072&label=PyPI&color=4B8BBE
    :target: https://pypi.org/project/ufjc/

.. |conda| image:: https://img.shields.io/conda/v/mrbuche/ufjc.svg?logo=anaconda&color=3EB049&label=Anaconda
    :target: https://anaconda.org/mrbuche/ufjc/

.. |docker| image:: https://img.shields.io/docker/v/mrbuche/ufjc?color=0db7ed&label=Docker%20Hub&logo=docker&logoColor=0db7ed
    :target: https://hub.docker.com/r/mrbuche/ufjc

.. |license| image:: https://img.shields.io/github/license/sandialabs/ufjc?label=License
    :target: https://github.com/sandialabs/ufjc/blob/main/LICENSE

.. |zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.6114263.svg
    :target: https://doi.org/10.5281/zenodo.6114263
