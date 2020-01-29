pyeemd – a Python interface to libeemd
======================================

This program is a Python interface to [libeemd][], which is a C library for
performing the ensemble empirical mode decomposition (EEMD), its complete
variant (CEEMDAN) or the regular empirical mode decomposition (EMD). The
details of what libeemd actually computes are available as a separate
[article][], which you should read if you are unsure about what EMD, EEMD and
CEEMDAN are.

[article]: https://dx.doi.org/10.1007/s00180-015-0603-9
[libeemd]: https://bitbucket.org/luukko/libeemd

Program license
---------------

pyeemd is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Installation
------------

### Dependencies

To use pyeemd you'll first need to download and install [libeemd][]. At runtime
pyeemd will try to locate the libeemd library (file `libeemd.so`) using
`ctypes.util.find_library` from the Python standard library. If you have
trouble getting pyeemd to find libeemd, you can also add `libeemd.so` to the
same directory as the file `pyeemd.py`.

To install pyeemd you'll also need
[setuptools](https://pypi.python.org/pypi/setuptools). Other dependencies are
handled automatically by setuptools.

pyeemd should work with either Python 2 (2.6 or later) or Python 3. If it
doesn't, please file an [issue](https://bitbucket.org/luukko/pyeemd/issues).

### Installation

Simply run

	./setup.py install

if installing system-wide, or

	./setup.py install --user

if installing only for a specific user.

Unit tests can be run with

	./setup.py test --quiet

To install for a specific Python version use that version to run `setup.py`, e.g.,

	python3 setup.py install --user

Additional documentation
------------------------

For code examples, please see the `examples` subdirectory.

For additional documentation, please look at the `doc` subdirectory (you'll
need `Sphinx` and `numpydoc` to build the documentation) or head straight to
[Read the Docs](http://pyeemd.readthedocs.org/).

Developing pyeemd
-----------------

Please report any issues using the Bitbucket issue tracker. If submitting pull
requests, please use `develop` as a target branch.
