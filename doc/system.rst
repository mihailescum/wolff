
*************************************
Constructing a System & Running Wolff
*************************************

The core of the library lies in the :cpp:class:`system` class and its member functions. Here, the state of your model is stored and cluster flip Monte Carlo can be carried out in various ways.

Note that the member objects and functions described here will change when compiled with certain compiler flags active, as described in :ref:`compile`.

.. cpp:class:: template\<class R_t, class X_t> system

Member Objects
==============

.. cpp:member:: v_t system::nv

   Stores the number of ordinary sites in the model.

.. cpp:member:: v_t system::ne

   Stores the number of ordinary bonds in the model.

.. cpp:member:: graph system::G

   Stores the graph describing the lattice, including the ghost site.

.. cpp:member:: double system::T

   Stores the temperature of the model.

.. cpp:member:: std::vector<X_t> system::s

   The :math:`i\text{th}` component stores the spin state of site :math:`i`.

.. cpp:member:: R_t system::s0

   Stores the state of the ghost site.

.. cpp:member:: std::function <double(const X_t&, const X_t&)> system::Z

   The function that returns the coupling between neighboring sites.

.. cpp:member:: std::function <double(const X_t&)> system::B

   The function that returns the coupling to the field.

Member Functions
================

.. cpp:function:: system::system(graph G, double T, std::function <double(const X_t&, const X_t&)> Z, std::function <double(const X_t&)> B)

   The constructor for systems. The arguments given are a graph :cpp:any:`G` *without* the ghost spin added, the temperature :cpp:any:`T`, and the coupling functions :cpp:any:`Z` and :cpp:any:`B`. The states of the spins and ghost site are initialized using the default constructors for :cpp:type:`X_t` and :cpp:type:`R_t`, respectively. :cpp:any:`nv` and :cpp:any:`ne` are taken directly from :cpp:any:`G`, after which the ghost site is added to :cpp:any:`G`.

.. cpp:function:: system::flip_cluster(v_t i0, const R_t& r, std::mt19937& rng, measurement<R_t, X_t>& A)

   Performs one Wolff cluster flip to the system. The cluster is seeded at vertex :cpp:any:`i0` and the spins added are transformed by :cpp:any:`r`. A random number generator :cpp:any:`rng` provides required random numbers during, and the relevant measurement hooks defined in the inherited class of :cpp:any:`A` are run during the cluster formation.

.. cpp:function:: system::run_wolff(N_t N, std::function <R_t(std::mt19937&, const system<R_t, X_t>&, v_t)> r_gen, measurement<R_t, X_t>& A, std::mt19937& rng)

   Performs :cpp:any:`N` Wolff cluster flips to the system. One must provide a function :cpp:func:`r_gen` that takes a random number generator, the system state, and the index of the seed site and returns a transformation for each flip. Also required are an object from a class which inherits :cpp:class:`measurement` to run measurements, and a random number generator :cpp:any:`rng`.

