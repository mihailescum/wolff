
************
Introduction
************

A library for running the Wolff algorithm on arbitrary systems in arbitrary
fields. A "spin state" and "spin symmetry transformation" type must be supplied
to define a system, along with a spin-spin coupling, a spin-field coupling,
and a generator of transformations of rank two. The library then supplies the
tools to run Wolff cluster-flip Monte Carlo on the resulting system, with
arbitrary measurements taken along the way.

A detailed description of the algorithm and its requirements can be found at https://arxiv.org/abs/1805.04019.

Getting Wolff
=============

This source for this library is available at https://git.kent-dobias.com/wolff/.

Installation
============

The only dependencies are a modern C++ compiler, cmake, and the standard libraries. With those at hand, the Wolff library and the provided examples can be built and installed using

.. code-block:: bash

   git clone https://git.kent-dobias.com/wolff/
   mkdir wolff/build
   cd wolff/build
   cmake ..
   make install

Custom install paths and compiler optimizations can be passed to cmake in the standard ways.

License & Use
=============

Wolff is licensed under the MIT license, a copy of which is included with the source code. The terms can be found here: https://spdx.org/licenses/MIT.html.

If you use or modify the library for use for academic purposes, please cite *Cluster representations and the Wolff algorithm in arbitrary external fields*, Jaron Kent-Dobias & James P Sethna, arXiv:1805.04019 [cond-mat.stat-mech].

