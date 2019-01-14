
.. _finite_states:

*************
Finite States
*************

One of the slower steps in running the algorithm is taking the exponent involved in computing the bond activation probabilities for each prospective bond. When the spins in your system have a finite number of possible states, the algorithm can be sped up considerably by precomputing the bond activation probabilities for every possible pair of spins. Once the appropriate things have been defined for your model, the compile definition :c:macro:`WOLFF_USE_FINITE_STATES` can be set to automate this process. The provided model headers :file:`wolff_models/ising.hpp` and :file:`wolff_models/potts.hpp` demonstrate the expected usage.

Required Definitions
====================

.. c:macro:: WOLFF_USE_FINITE_STATES

   This macro must defined before :file:`wolff.hpp` or any of the other header files are invoked.

.. c:macro:: WOLFF_FINITE_STATES_N

   This macro must be defined and given a value equal to the number of states your model can take. It must be defined before :file:`wolff.hpp` or any of the other header files are invoked.

.. cpp:function:: X_t::X_t(q_t)

   Your spin class :cpp:class:`X_t` must have a constructor defined that takes a :cpp:type:`q_t` and returns a unique state for all arguments less than :c:macro:`WOLFF_FINITE_STATES_N`.

.. cpp:function:: q_t X_t::enumerate()

   Your spin class :cpp:class:`X_t` must have a function defined that returns the index associated with a given state. This must be the inverse function of the constructor above.

