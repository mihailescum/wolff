
.. default-domain:: cpp

************
Measurements
************

One of the most complicated parts of using this library is setting up measurements to take during a run, but once understood they provide a powerful way to quickly measure arbitrary properties. The class is defined in the header :file:`lib/wolff.hpp`.

.. class:: template<R_t, X_t, G_t = graph<>> measurement

   This class defines several virtual member functions. To use the library, you must supply your own measurement class that inherits this one and defines those functions, which may be trivial.

   Each member function defines a hook that is run at a specific time during execution. These hooks can be used to modify member objects of your inheritor measurement class, and thereby extract information from the simulation.

   .. function:: void measurement::pre_cluster(unsigned n, unsigned N, const system<R_t, X_t, G_t>& S, const G_t::vertex& v, const R_t& r)

      This hook is run prior to cluster formation.

      :param unsigned n: The number of clusters already flipped.
      :param unsigned N: The total number of clusters to flip.
      :param const system<R_t, X_t, G_t>& S: The current system state.
      :param const G_t\:\:vertex& v: The seed vertex for the cluster.
      :param const R_t& r: The transformation by which the cluster will be flipped.

   .. function:: void measurement::plain_bond_visited(const system<R_t, X_t, G_t>& S, const G_t::halfedge& e, const X_t& s1_new, double dE)

         This hook is run when an ordinary edge is visited during cluster formation, whether it is activated or not.

         :param const system<R_t, X_t, G_t>& S: The current system state.
         :param const G_t\:\:halfedge& e:   The edge being considered, with e.self referencing the vertex soon to be flipped.
         :param const X_t& s1_new:         The state the vertex will be flipped to.
         :param double dE:                 The change in energy that will result when the spin is flipped.

   .. cpp:function:: void measurement::plain_site_transformed(const system<R_t, X_t, G_t>& S, const G_t::vertex& v, const X_t& si_new)

         This hook is run when an ordinary site is about to be flipped.

         :param const system<R_t, X_t, G_t>& S: The current system state.
         :param const G_t\:\:vertex& v:                     The vertex soon to be flipped.
         :param const X_t& si_new:          The state that vertex will be flipped to.

   .. cpp:function:: void measurement::ghost_bond_visited(const system<R_t, X_t, G_t>& S, const G_t::vertex& v, const X_t& s0si_old, const X_t& s0si_new, double dE)

         This hook is run when an edge containing the ghost site is visited during cluster formation, whether activated or not.

         :param const system<R_t, X_t, G_t>& S: The current system state.
         :param const G_t\:\:vertex& v:                     The non-ghost site in this edge.
         :param const X_t& s0si_old:       The state that the spin on the non-ghost site has before the transformation is applied, rotated by the inverse action of the ghost site.
         :param const X_t& s0si_new:       The state that the spin on the non-ghost site will have after the transformation is applied, rotated by the inverse action of the ghost site.
         :param double dE:                 The change in energy that will result when one site is transformed.

   .. cpp:function:: void measurement::ghost_site_transformed(const system<R_t, X_t, G_t>& S, const R_t& R_new)

         This hook is run when the ghost site is about to be flipped.

         :param const system<R_t, X_t, G_t>& S: The current system state.
         :param const R_t& R_new:          The state the ghost site will be flipped to.

   .. cpp:function:: void measurement::post_cluster(unsigned n, unsigned N, const system<R_t, X_t, G_t>& S) 

         This hook is run after a cluster has been flipped.

         :param unsigned n:                     The number of clusters already flipped.
         :param unsigned N:                     The total number of clusters to flip.
         :param const system<R_t, X_t, G_t>& S: The current system state.

