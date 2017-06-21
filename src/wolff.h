
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <float.h>
#include <sys/types.h>
#include <inttypes.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <stdbool.h>
#include <assert.h>
#include <fftw3.h>

#include <jst/graph.h>
#include <jst/rand.h>

typedef struct {
  graph_t *g;
  bool *spins;
  int32_t M;
  double H;
} ising_state_t;

typedef struct ll_tag {
	uint32_t x;
	struct ll_tag *next;
} ll_t;

typedef struct {
  uint32_t nv;
  double dH;
  int32_t dM;
} cluster_t;

double get_hamiltonian(graph_t *g, double *coupling, bool *x);

void stack_push(ll_t **q, uint32_t x);

uint32_t stack_pop(ll_t **q);

bool stack_contains(const ll_t *q, uint32_t x);

cluster_t *flip_cluster(const graph_t *g, const double *ps, double H, bool *x, gsl_rng *r);

graph_t *graph_add_ext(const graph_t *g);

double hh(double th);

double *get_bond_probs(double T, double H, ising_state_t *s);

uint32_t wolff_step(double T, double H, ising_state_t *s, gsl_rng *r, double *ps);

