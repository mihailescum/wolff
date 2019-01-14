
.. default-domain:: cpp

*************************************
Constructing a System & Running Wolff
*************************************

This class and associated member functions are defined in the header file :file:`lib/wolff.hpp`.

.. class:: template\<R_t, X_t, G_t = graph\<>> system

   The core of the library lies in the :class:`system` class and its member functions. Here, the state of your model is stored and cluster flip Monte Carlo can be carried out in various ways.

   Note that the member objects and functions described here will change when compiled with certain compiler flags active, as described in :ref:`compile`.

   .. member:: unsigned system::nv

      Stores the number of ordinary sites in the model.

   .. member:: unsigned system::ne

      Stores the number of ordinary bonds in the model.

   .. member:: G_t system::G

      Stores the graph describing the lattice, including the ghost site. See :ref:`graphs` for details on this class. The template parameter has a default option corresponding to a graph type with no special vertex or edge properties.

   .. member:: double system::T

      Stores the temperature of the model.

   .. member:: std::vector<X_t> system::s

      The :math:`i\text{th}` component stores the spin state of site :math:`i`.

   .. member:: R_t system::s0

      Stores the state of the ghost site.

   .. member:: std::function <double(const X_t&, const X_t&)> system::Z

      The function that returns the coupling between neighboring sites.

   .. member:: std::function <double(const X_t&)> system::B

      The function that returns the coupling to the field.

   .. function:: system::system(G_t G, double T, std::function <double(const X_t&, const X_t&)> Z, std::function <double(const X_t&)> B)

         The constructor for systems.
         
         :param G_t G: A lattice graph *without* the ghost spin added.
         :param double T: The temperature.
         :param std\:\:function<double(const X_t&, const X_t&)> Z: The bond coupling.
         :param std\:\:function<double(const X_t&)> B: The field coupling.
         
         The states of the spins and ghost site are initialized using the default constructors for :type:`X_t` and :type:`R_t`, respectively. :any:`nv` and :any:`ne` are taken directly from :any:`G`, after which the ghost site is added to :any:`G`.

   .. function:: system::flip_cluster(const G_t::vertex& v0, const R_t& r, std::mt19937& rng, measurement<R_t, X_t, G_t>& A)

         Performs one Wolff cluster flip to the system. 
         
         :param const G_t\:\:vertex& v0: The vertex of the seed site.
         :param const R_t& r: The transformation by which the cluster is flipped.
         :param std\:\:mt19937& rng: A random number generator.
         :param measurement<R_t, X_t, G_t>& A: Object whose class inherits :class:`measurement` and provides relevant measurement hooks.

   .. function:: system::run_wolff(unsigned N, std::function <R_t(std::mt19937&, const system<R_t, X_t, G_t>&, const G_t::vertex&)> r_gen, measurement<R_t, X_t, G_t>& A, std::mt19937& rng)

         Performs :any:`N` Wolff cluster flips to the system.
         
         :param unsigned N: Number of clusters to flip.
         :param std\:\:function <R_t(std\:\:mt19937&, const system<R_t, X_t, G_t>&, const G_t\:\:vertex&>)> r_gen: Generator of transformations for the cluster flips.
         :param measurement<R_t, X_t, G_t>& A: Object whose class inherits :class:`measurement` and provides relevant measurement hooks.
         :param std\:\:mt19937& rng: A random number generator.


