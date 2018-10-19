
************
Measurements
************

One of the most complicated parts of using this library is setting up measurements to take during a run, but once understood they provide a powerful way to quickly measure arbitrary properties. The class is defined in the header :file:`wolff/measurement.hpp`.

.. cpp:class:: template\<class R_t, class X_t> measurement

   This class defines several virtual member functions. To use the library, you must supply your own measurement class that inherits this one and defines those functions, which may be trivial.

   Each member function defines a hook that is run at a specific time during execution. These hooks can be used to modify member objects of your inheritor measurement class, and thereby extract information from the simulation.

.. cpp:function:: template<void measurement::pre_cluster(N_t n, N_t N, const system<R_t, X_t>& S, v_t i0, const R_t& r)

   This hook is run prior to cluster formation. The arguments supplied are:

   =========================== ===
   :cpp:var:`n`                the number of clusters already flipped
   :cpp:var:`N`                the total number of clusters to flip
   :cpp:var:`S`                the current system state
   :cpp:var:`i0`               the index of the seed vertex for the cluster
   :cpp:var:`r`                the transformation by which the cluster will be flipped
   =========================== ===

.. cpp:function:: void measurement::plain_bond_visited(const system<R_t, X_t>& S, v_t i1, const X_t& s1_new, v_t i2, double dE)

   This hook is run when an ordinary edge is visited during cluster formation, whether it is activated or not. The arguments supplied are:

   ================================= ===
   :code:`const system<R_t, X_t>& S` the current system state
   :code:`v_t i1`                    the index of the vertex soon to be flipped
   :code:`const X_t& s1_new`         the state vertex :cpp:var:`i1` will be flipped to
   :code:`v_t i2`                    the other vertex sharing the edge
   :code:`double dE`                 the change in energy that will result when the spin at site :cpp:var:`i1` is flipped
   ================================= ===

.. cpp:function:: void measurement::plain_site_transformed(const system<R_t, X_t>& S, v_t i, const X_t& si_new)

   This hook is run when an ordinary site is about to be flipped. The arguments supplied are:

   ================================= ===
   :code:`const system<R_t, X_t>& S` the current system state
   :code:`v_t i`                     the index of the vertex soon to be flipped
   :code:`const X_t& si_new`         the state vertex :cpp:var:`i` will be flipped to
   ================================= ===

|

.. cpp:function:: void measurement::ghost_bond_visited(const system<R_t, X_t>& S, v_t i, const X_t& s0si_old, const X_t& s0si_new, double dE)

   This hook is run when an edge containing the ghost site is visited during cluster formation, whether activated or not. The arguments supplied are:

   ================================= ===
   :code:`const system<R_t, X_t>& S` the current system state
   :code:`v_t i`                     the index of the non-ghost site in this edge
   :code:`const X_t& s0si_old`       the state that the spin on the non-ghost site has before the transformation is applied, rotated by the inverse action of the ghost site
   :code:`const X_t& s0si_new`       the state that the spin on the non-ghost site will have after the transformation is applied, rotated by the inverse action of the ghost site
   :code:`double dE`                 the change in energy that will result when one site is transformed
   ================================= ===

|

.. cpp:function:: void measurement::ghost_site_transformed(const system<R_t, X_t>& S, const R_t& R_new)

   This hook is run when the ghost site is about to be flipped. The arguments supplied are:

   ================================= ===
   :code:`const system<R_t, X_t>& S` the current system state
   :code:`const R_t& R_new`          the state the ghost site will be flipped to
   ================================= ===

|

.. cpp:function:: void measurement::post_cluster(N_t n, N_t N, const system<R_t, X_t>& S) 

   This hook is run after a cluster has been flipped. The arguments supplied are:

   ================================= ===
   :code:`N_t n`                     the number of clusters already flipped
   :code:`N_t N`                     the total number of clusters to flip
   :code:`const system<R_t, X_t>& S` the current system state
   ================================= ===

