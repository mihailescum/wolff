
************
Custom Types
************

The Wolff library uses several custom types, defined to promote memory and cache effiency, promote shorter lines, and associate natural scope to certain parameters. They are all aliases for `standard C99 fixed width integer types`_. These are defined in the header file :file:`wolff/types.h`. Here is a list of the types and what they are used to hold:

.. cpp:type:: uint_fast32_t v_t

   Holds indicies for lattice vertices and edges.

.. cpp:type:: uint_fast8_t q_t

   Holds indicies for enumerating states, as in :math:`q`-state Potts.

.. cpp:type:: uint_fast8_t D_t

   Holds the dimension of space.

.. cpp:type:: uint_fast16_t L_t

   Holds the linear extent of the lattice.

.. cpp:type:: uint_fast64_t N_t

   Holds a count of Monte Carlo steps.

Other Definitions
=================

All types are unsigned integers of a minimum size appropriate for most realistic purposes, with the :c:type:`fast` directives meaning that your compiler may choose a larger type for efficiency's sake.

For each type :c:type:`x_t`, the macro :c:macro:`MAX_x` supplies the maximum value the type may hold, the macro :c:macro:`PRIx` supplies the format constants for use with the :c:func:`fprintf` family of functions, and the macro :c:macro:`SCNx` supplies the format constants for use with the :c:func:`fscanf` family of functions.

.. _standard C99 fixed width integer types: https://en.cppreference.com/w/c/types/integer

