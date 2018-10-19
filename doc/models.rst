
****************
Defining a Model
****************

To define your own model, you need a class for spin states, a class for spin transformations, a bond and site coupling, and a generator of transformations. Many examples of models can be found in the headers contained in :file:`wolff/models/`.
 
Spin States
===========

.. cpp:class:: X_t

   The spin type, which we have been denoting with :cpp:class:`X_t`, only needs to have a default constructor. If your spins can take only finitely many values, consider following the instructions in :ref:`finite_states` to significantly speed the algorithm.

Transformations
===============

.. cpp:class:: R_t

   The transformation type must have a default constructor, and the following member functions:

   .. cpp:function:: X_t act(const X_t&)

      The action of the transformation on a spin, returning the transformed spin.

      :param const X_t&: The spin state to transform.
      :returns: The transformed spin state.

   .. cpp:function:: R_t act(const R_t&)

      The action of the transformation on another transformation, returning the transformed transformation.

      :param const R_t&: The transformation state to transform.
      :returns: The transformed spin state.

   .. cpp:function:: X_t act_inverse(const X_t&)

      The action of the inverse transformation on a spin, returning the transformed spin.

      :param const X_t&: The spin state to transform.
      :returns: The transformed spin state.

   .. cpp:function:: R_t act_inverse(const R_t&)

      The action of the inverse transformation on another transformation, returning the transformed transformation.

      :param const R_t&: The transformation state to transform.
      :returns: The transformed spin state.

Couplings
=========

When a :cpp:class:`system` object is initialized it needs to be given a bond and field coupling, to resemble the Hamiltonian

.. math::
         \mathcal H = -\sum_{\langle ij\rangle}Z(s_i,s_j)-\sum_iB(s_i)

Note that building with certain compile flags will change the form that these coupling functions must take, as outlined in :ref:`compile`.

Bond Coupling
-------------

.. cpp:function:: double Z(const X_t&, const X_t&)

   A function that takes two spins and outputs the coupling between them. For traditional spin models, this is typically something like a dot product.

Field Coupling
--------------

.. cpp:function:: double B(const X_t&)

   A function that takes one spin and outputs the coupling between it and an external field.

