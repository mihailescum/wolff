
.. _finite_states:

*************
Finite States
*************

One of the slower steps in running the algorithm is taking the exponent involved in computing the bond activation probabilities for each prospective bond. When the spins in your system have a finite number of possible states, the algorithm can be sped up considerably by precomputing the bond activation probabilities for every possible pair of spins. Once the appropriate things have been defined for your model, the header :file:`wolff/finite_states.hpp` can be invoked to automate this process. The provided model headers :file:`wolff/models/ising.hpp` and :file:`wolff/models/potts.hpp` demonstrate the expected usage.

Required Definitions
====================

.. c:macro:: WOLFF_FINITE_STATES_N

   This macro must be defined and given a value equal to the number of states your model can take.

.. cpp:var:: const X_t finite_states_possible[WOLFF_FINITE_STATES_N]

   You must supply a constant array which lists each of the possible states of your individual spins. 

.. cpp:function:: q_t finite_states_enum(const X_t&)

   You must define a function which takes each state and returns its index in :cpp:var:`finite_states_possible`.

