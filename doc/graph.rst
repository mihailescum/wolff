
.. default-domain:: cpp

******
Graphs
******

This class is defined in the header :file:`wolff/graph.hpp`.

.. class:: graph

   Lattices are described by objects of class :class:`graph`, a minimal implementation of graphs.

   .. member:: D_t graph::D

      The dimension of the graph. This property is unused by the core library.

   .. member:: L_t graph::L

      The linear size of the graph. This property is unused by the core library.

   .. member:: v_t graph::ne

      The number of edges in the graph. This property is unused by the core library.

   .. member:: v_t graph::nv

      The number of vertices in the graph.

   .. member:: std::vector<std::vector<v_t>> graph::adjacency

      The adjacency list for the graph. The :math:`i\text{th}` element of :member:`adjacency` is a standard library vector containing the indices of all vertices adjacent to vertex :math:`i`.

   .. member:: std::vector<std::vector<double>> graph::coordinate

      The coordinates of the graph vertices. The :math:`i\text{th}` element of :var:`coordinate` is a standard library vector of length :var:`D` containing the spatial coordinates of vertex :math:`i`. This property is unused by the core library.

   .. function:: graph::graph()

      The default constructor. Initializes an empty graph, i.e., :var:`D`, :var:`L`, :var:`ne`, and :var:`nv` are all zero and :var:`adjacency` and :var:`coordinate` are uninitialized.

   .. function:: graph::graph(D_t D, L_t L)

      Initializes a graph of a :var:`D`-dimensional hypercubic lattice with :var:`L` vertices per side. This is the only nontrivial graph constructor supplied by the core library. The library will work with arbitrary graphs, and if a different lattice is needed consider calling the default constructor and populating the member objects youself before handing the graph to the :class:`system` constructor.

      :param D_t D: The dimension of space.
      :param L_t L: The number of vertices per edge of the hypercube.

   .. function:: void graph::add_ghost()

      Calling this function on a graph object will add a ghost site to the graph. Explicitly, a new vertex is added that is adjacent to every other vertex in the graph. This vertex will have the last index, which is equal to number of vertices in the original graph.

