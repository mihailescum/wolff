
.. _compile:

***************
Compile Options
***************

There are several features that are not active by default, but can be used by defining a compiler flag. These will modify the constructor of :cpp:class:`system`. These can be used in concert when applicable. Note that at the moment, there is no :ref:`finite_states` support for systems with bond or site dependence.

Building Without A Field
========================

.. c:macro:: WOLFF_NO_FIELD

   When this flag is defined Wolff will not use a field. When a field isn't necessary, this will improve performance: no ghost site will be initialize and no time will be wasted checking the energy change with respect to a uncoupled ghost site. 

   When defined, the constructor for :cpp:class:`system` becomes :cpp:any:`system(graph, double, std::function <double(const X_t&, const X_t&)>)`. The resulting :cpp:class:`system` object does not have member objects :cpp:any:`B` or :cpp:member:`s0`, and :cpp:member:`system::G` does not have a ghost site initialized.

Building With Bond Dependence
=============================

.. c:macro:: WOLFF_BOND_DEPENDENCE

   When this flag is defined Wolff will ask for the indices of the spins on each side a bond, allowing the implementation of random bonds or anisotropic interactions.

   When defined, the bond coupling must be a function of the form :cpp:any:`double Z(const G_t::halfedge&, const X_t&, const X_t&)`, where the first argument is the edge being considered. A function of this type is passed to :cpp:class:`system` in place of the original bond coupling.


Building With Site Dependence
=============================

.. c:macro:: WOLFF_SITE_DEPENDENCE

   When this flag is defined Wolff will ask for the indices of the spin when measuring the external field, allowing the implementation of random fields or to emulate boundaries.

   When defined, the field coupling must be a function of the form :cpp:any:`double B(const G_t::vertex&, const X_t&)`, where the first argument is the vertex the spin is on. A function of this type is passed to :cpp:class:`system` in place of the original field coupling.

   An example of a system of this type can be found in :file:`examples/ising_random_field.cpp`, which uses a non-trivial vertex property to communicate vertex dependence to the field coupling function.
