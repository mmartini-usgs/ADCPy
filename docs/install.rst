Installing ADCPy
****************

If you do not already have python installed, we recommend Anaconda.  In particular we use this installation:  http://ioos.github.io/notebooks_demos/other_resources/

Using git, clone ADCPy from github at:  https://github.com/mmartini-usgs/ADCPy.

``cd`` to the directory containing the package, and type:

``pip install -e . --no-deps``

This will not install dependencies. If you get errors about missing packages, install them using ``conda install`` or ``pip install``.

Note also that you will import one of the three sub-modules as follows:

* ``import EPICstuff``
* ``import TRDIstuff``
* ``import Nortekstuff``
