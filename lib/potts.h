#pragma once

#include <cmath>
#include <stdio.h>

#include "types.h"
#include "vector.h"

/* The following is the minimum definition of a spin class.
 *
 * The class must contain an M_t and an F_t for holding the sum of an
 * integer number of spins and a double-weighted number of spins,
 * respectively.
 *
 * void init(X_t *p);
 * void free_spin(X_t p);
 * void free_spin(M_t p);
 * void free_spin(F_t p);
 * X_t copy(X_t x);
 * void add(M_t *x1, int factor, X_t x2);
 * void add(F_t *x1, double factor, X_t x2);
 * M_t scalar_multiple(int factor, X_t x);
 * F_t scalar_multiple(double factor, X_t x);
 * double norm_squared(F_t x);
 * void write_magnetization(M_t M, FILE *outfile);
 *
 */

template <q_t q>
class potts_t {
  public:
    q_t x;

    typedef vector_t<q, int> M_t;
    typedef vector_t<q, double> F_t;
};

template <q_t q>
void init(potts_t <q> *p) {
  p->x = 0;
}

template <q_t q>
void free_spin(potts_t <q> s) {
  // do nothing!
}

template <q_t q>
potts_t <q> copy(potts_t <q> s) {
  return s;
}

template <q_t q>
void add(vector_t<q, int> *s1, int a, potts_t <q> s2) {
  (s1->x)[s2.x] += a;
}

template <q_t q>
void add(vector_t<q, double> *s1, double a, potts_t <q> s2) {
  (s1->x)[s2.x] += a;
}

template <q_t q>
vector_t<q, int> scalar_multiple(int factor, potts_t <q> s) {
  vector_t<q, int> M;
  M.x = (int *)calloc(q, sizeof(int));
  M.x[s.x] += factor;
  return M;
}

template <q_t q>
vector_t<q, double> scalar_multiple(double factor, potts_t <q> s) {
  vector_t<q, double> F;
  F.x = (double *)calloc(q, sizeof(double));
  F.x[s.x] += factor;
  return F;
}

// we could inherit norm_squared from vector.h, but convention dictates that
// potts norms be changed by a constant factor
template <q_t q>
double norm_squared(vector_t<q, double> s) {
  double total = 0;
  for (q_t i = 0; i < q; i++) {
    total += pow(s.x[i], 2);
  }

  return total * (double)q / ((double)q - 1.0);
}

// we could inherit write_magnetization from vector.h, but since M.x must sum
// to nv we don't need to write the last element
template <q_t q>
void write_magnetization(vector_t<q, int> M, FILE *outfile) {
  fwrite(M.x, sizeof(int), q - 1, outfile);
}

// knock yourself out
const potts_t<POTTSQ> states[256] = {{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {11}, {12}, {13}, {14}, {15}, {16}, {17}, {18}, {19}, {20}, {21}, {22}, {23}, {24}, {25}, {26}, {27}, {28}, {29}, {30}, {31}, {32}, {33}, {34}, {35}, {36}, {37}, {38}, {39}, {40}, {41}, {42}, {43}, {44}, {45}, {46}, {47}, {48}, {49}, {50}, {51}, {52}, {53}, {54}, {55}, {56}, {57}, {58}, {59}, {60}, {61}, {62}, {63}, {64}, {65}, {66}, {67}, {68}, {69}, {70}, {71}, {72}, {73}, {74}, {75}, {76}, {77}, {78}, {79}, {80}, {81}, {82}, {83}, {84}, {85}, {86}, {87}, {88}, {89}, {90}, {91}, {92}, {93}, {94}, {95}, {96}, {97}, {98}, {99}, {100}, {101}, {102}, {103}, {104}, {105}, {106}, {107}, {108}, {109}, {110}, {111}, {112}, {113}, {114}, {115}, {116}, {117}, {118}, {119}, {120}, {121}, {122}, {123}, {124}, {125}, {126}, {127}, {128}, {129}, {130}, {131}, {132}, {133}, {134}, {135}, {136}, {137}, {138}, {139}, {140}, {141}, {142}, {143}, {144}, {145}, {146}, {147}, {148}, {149}, {150}, {151}, {152}, {153}, {154}, {155}, {156}, {157}, {158}, {159}, {160}, {161}, {162}, {163}, {164}, {165}, {166}, {167}, {168}, {169}, {170}, {171}, {172}, {173}, {174}, {175}, {176}, {177}, {178}, {179}, {180}, {181}, {182}, {183}, {184}, {185}, {186}, {187}, {188}, {189}, {190}, {191}, {192}, {193}, {194}, {195}, {196}, {197}, {198}, {199}, {200}, {201}, {202}, {203}, {204}, {205}, {206}, {207}, {208}, {209}, {210}, {211}, {212}, {213}, {214}, {215}, {216}, {217}, {218}, {219}, {220}, {221}, {222}, {223}, {224}, {225}, {226}, {227}, {228}, {229}, {230}, {231}, {232}, {233}, {234}, {235}, {236}, {237}, {238}, {239}, {240}, {241}, {242}, {243}, {244}, {245}, {246}, {247}, {248}, {249}, {250}, {251}, {252}, {253}, {254}, {255}};
template <q_t q>
q_t state_to_ind(potts_t<q> state) { return (q_t)state.x; }

