
.. default-domain:: cpp
.. _graphs:

******
Graphs
******

This class is defined in the header :file:`lib/wolff.hpp`.

.. class:: \template <class vertex_prop = std::tuple\<>, class edge_prop = std::tuple\<>> graph

   Lattices are described by objects of class :class:`graph`, a minimal implementation of graphs. Can be called with `graph<>` if no properties need to be associated with vertices or edges. Otherwise, those properties can be supplied as classes.

   .. member:: unsigned graph::D

      The dimension of the graph. This property is unused by the core library.

   .. member:: unsigned graph::L

      The linear size of the graph. This property is unused by the core library.

   .. member:: unsigned graph::ne

      The number of edges in the graph. This property is unused by the core library.

   .. member:: unsigned graph::nv

      The number of vertices in the graph.

   .. class:: graph::vertex

      This class describes the vertices on the graph.

      .. member:: unsigned ind

         The index of the vertex, which is also its position in the list :var:`vertices`.

      .. member:: std::list<halfedge> edges

         The list of edges emanating from this vertex.

      .. member:: vertex_prop prop

         Template-defined class which stores optional properties of the vertex.

   .. class:: graph::halfedge

      This class describes the halfedges on the graph.

      .. member:: vertex& self

         A reference to the vertex this halfedge is emanating from.

      .. member:: vertex& neighbor

         A reference to the vertex this halfedge is going to.

      .. member:: edge_prop prop

         Template-defined class which stores optional properties of the edge.

      .. function:: halfedge(vertex &self, vertex &neighbor)

         Constructor which sets self and neighbor from supplied vertices.

   .. member:: std::vector<vertex> graph::vertices

      A list of all vertices in the graph.

   .. function:: graph::graph()

      The default constructor. Initializes an empty graph, i.e., :var:`D`, :var:`L`, :var:`ne`, and :var:`nv` are all zero and :var:`vertices` is uninitialized.

   .. function:: graph::graph(unsigned D, unsigned L)

      Initializes a graph of a :var:`D`-dimensional hypercubic lattice with :var:`L` vertices per side. This is the only nontrivial graph constructor supplied by the core library. The library will work with arbitrary graphs, and if a different lattice is needed consider calling the default constructor and populating the member objects youself before handing the graph to the :class:`system` constructor.

      :param unsigned D: The dimension of space.
      :param unsigned L: The number of vertices per edge of the hypercube.

   .. function:: void graph::add_ghost()

      Calling this function on a graph object will add a ghost site to the graph. Explicitly, a new vertex is added that is adjacent to every other vertex in the graph. This vertex will have the last index, which is equal to number of vertices in the original graph.

