
******
Graphs
******

This class is defined in the header :file:`wolff/graph.hpp`.

.. cpp:class:: graph

   Lattices are described by objects of class :cpp:class:`graph`, a minimal implementation of graphs.

.. cpp:member:: D_t graph::D

   The dimension of the graph. This property is unused by the core library.

.. cpp:member:: L_t graph::L

   The linear size of the graph. This property is unused by the core library.

.. cpp:member:: v_t graph::ne

   The number of edges in the graph. This property is unused by the core library.

.. cpp:member:: v_t graph::nv

   The number of vertices in the graph.

.. cpp:member:: std::vector<std::vector<v_t>> graph::adjacency

   The adjacency list for the graph. The :math:`i\text{th}` element of :cpp:member:`adjacency` is a standard library vector containing the indices of all vertices adjacent to vertex :math:`i`.

.. cpp:member:: std::vector<std::vector<double>> graph::coordinate

   The coordinates of the graph vertices. The :math:`i\text{th}` element of :cpp:var:`coordinate` is a standard library vector of length :cpp:var:`D` containing the spatial coordinates of vertex :math:`i`. This property is unused by the core library.

.. cpp:function:: graph::graph()

   The default constructor. Initializes an empty graph, i.e., :cpp:var:`D`, :cpp:var:`L`, :cpp:var:`ne`, and :cpp:var:`nv` are all zero and :cpp:var:`adjacency` and :cpp:var:`coordinate` are uninitialized.

.. cpp:function:: graph::graph(D_t D, L_t L)

   Calling the constructor with :cpp:var:`D` and :cpp:var:`L` will initialize the graph of a :cpp:var:`D`-dimensional hypercubic lattice with :cpp:var:`L` vertices per side. This is the only nontrivial graph constructor supplied by the core library.

.. cpp:function:: void graph::add_ghost()

   Calling this function on a graph object will add a ghost site to the graph. Explicitly, a new vertex is added that is adjacent to every other vertex in the graph. This vertex will have the last index, which is equal to number of vertices in the original graph.

