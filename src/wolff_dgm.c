
#include <getopt.h>

#include <cluster.h>

double identity(h_t x) {
  return -pow(x, 2);
}

double basic_H(double *H, h_t x) {
  return -H[0] * pow(x, 2);
}

int main(int argc, char *argv[]) {

  L_t L = 128;
  count_t N = (count_t)1e7;
  count_t min_runs = 10;
  count_t n = 3;
  D_t D = 2;
  double T = 2.26918531421;
  double *H = (double *)calloc(MAX_Q, sizeof(double));
  double eps = 0;
  bool silent = false;
  bool record_autocorrelation = false;
  count_t ac_skip = 1;
  count_t W = 10;

  int opt;
  q_t H_ind = 0;

  while ((opt = getopt(argc, argv, "N:n:D:L:T:H:m:e:saS:W:")) != -1) {
    switch (opt) {
    case 'N':
      N = (count_t)atof(optarg);
      break;
    case 'n':
      n = (count_t)atof(optarg);
      break;
    case 'D':
      D = atoi(optarg);
      break;
    case 'L':
      L = atoi(optarg);
      break;
    case 'T':
      T = atof(optarg);
      break;
    case 'H':
      H[H_ind] = atof(optarg);
      H_ind++;
      break;
    case 'm':
      min_runs = atoi(optarg);
      break;
    case 'e':
      eps = atof(optarg);
      break;
    case 's':
      silent = true;
      break;
    case 'a':
      record_autocorrelation = true;
      break;
    case 'S':
      ac_skip = (count_t)atof(optarg);
      break;
    case 'W':
      W = (count_t)atof(optarg);
      break;
    default:
      exit(EXIT_FAILURE);
    }
  }

  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, rand_seed());

  dgm_state_t *s = (dgm_state_t *)calloc(1, sizeof(dgm_state_t));

  graph_t *h = graph_create_square(D, L);
  s->g = graph_add_ext(h);

  s->spins = (h_t *)calloc(h->nv, sizeof(h_t));

  s->H_info = H;
  s->T = T;
  s->H = basic_H;
  s->J = identity;

  s->R = (dihinf_t *)calloc(1, sizeof(dihinf_t));

  s->M = 0;
  s->E = 0;

  double diff = 1e31;
  count_t n_runs = 0;
  count_t n_steps = 0;

  meas_t *E, *clust, *M, *dM;

  M = (meas_t *)calloc(1, sizeof(meas_t ));
  dM = (meas_t *)calloc(1, sizeof(meas_t ));

  E = calloc(1, sizeof(meas_t));
  clust = calloc(1, sizeof(meas_t));

  autocorr_t *autocorr;
  if (record_autocorrelation) {
    autocorr = (autocorr_t *)calloc(1, sizeof(autocorr_t));
    autocorr->W = 2 * W + 1;
    autocorr->OO = (double *)calloc(2 * W + 1, sizeof(double));
  }

  if (!silent) printf("\n");
  while (((diff > eps || diff != diff) && n_runs < N) || n_runs < min_runs) {
    if (!silent) printf("\033[F\033[JWOLFF: sweep %" PRIu64
           ", dH/H = %.4f, dM/M = %.4f, dC/C = %.4f, dX/X = %.4f, cps: %.1f\n",
           n_runs, fabs(meas_dx(E) / E->x), meas_dx(M) / M->x, meas_dc(E) / meas_c(E), meas_dc(M) / meas_c(M), h->nv / clust->x);

    count_t n_flips = 0;

    while (n_flips / h->nv < n) {
      v_t v0 = gsl_rng_uniform_int(r, h->nv);
      h_t step = round((((double)s->M) / h->nv) + gsl_ran_gaussian(r, 5));

      v_t tmp_flips = flip_cluster_dgm(s, v0, step, r);
      n_flips += tmp_flips;

      if (n_runs > 0) {
        n_steps++;
        meas_update(clust, tmp_flips);
      }

      if (record_autocorrelation && n_runs > 0) {
        if (n_steps % ac_skip == 0) {
          update_autocorr(autocorr, s->E);
        }
      }
    }

    meas_update(M, s->M);
    h_t min_h, max_h;
    min_h = MAX_H;
    max_h = MIN_H;
    for (v_t i = 0; i < h->nv; i++) {
      if (s->spins[i] < min_h) {
        min_h = s->spins[i];
      } else if (s->spins[i] > max_h) {
        max_h = s->spins[i];
      }
    }
    meas_update(dM, max_h - min_h);
    meas_update(E, s->E);

    diff = fabs(meas_dc(E) / meas_c(E));

    n_runs++;
  }
  if (!silent) {
    printf("\033[F\033[J");
  }
  printf("WOLFF: sweep %" PRIu64
         ", dH/H = %.4f, dM/M = %.4f, dC/C = %.4f, dX/X = %.4f, cps: %.1f\n",
         n_runs, fabs(meas_dx(E) / E->x), meas_dx(M) / M->x, meas_dc(E) / meas_c(E), meas_dc(M) / meas_c(M), h->nv / clust->x);

  double tau = 0;
  bool tau_failed = false;

  if (record_autocorrelation) {
    double *Gammas = (double *)malloc((W + 1) * sizeof(double));

    Gammas[0] = 1 + rho(autocorr, 0);
    for (uint64_t i = 0; i < W; i++) {
      Gammas[1 + i] = rho(autocorr, 2 * i + 1) + rho(autocorr, 2 * i + 2);
    } 

    uint64_t n;
    for (n = 0; n < W + 1; n++) {
      if (Gammas[n] <= 0) {
        break;
      }
    }

    if (n == W + 1) {
      printf("WARNING: correlation function never hit the noise floor.\n");
      tau_failed = true;
    }

    if (n < 2) {
      printf("WARNING: correlation function only has one nonnegative term.\n");
      tau_failed = true;
    }

    double *conv_Gamma = get_convex_minorant(n, Gammas);

    double ttau = - 0.5;

    for (uint64_t i = 0; i < n + 1; i++) {
      ttau += conv_Gamma[i];
    }
    
    free(Gammas);
    free(autocorr->OO);
    while (autocorr->Op != NULL) {
      stack_pop_d(&(autocorr->Op));
    }
    free(autocorr);
    
    tau = ttau * ac_skip * clust->x / h->nv;
  }

  if (tau_failed) {
    tau = 0;
  }

  FILE *outfile = fopen("out.m", "a");

  fprintf(outfile, "<|D->%" PRID ",L->%" PRIL ",T->%.15f", D, L, T);
  fprintf(outfile, ",E->%.15f,\\[Delta]E->%.15f,C->%.15f,\\[Delta]C->%.15f,M->%.15f", E->x / h->nv, meas_dx(E) / h->nv, meas_c(E) / h->nv, meas_dc(E) / h->nv, M->x / h->nv);
  fprintf(outfile, ",\\[Delta]M->%.15f", meas_dx(M) / h->nv);
  fprintf(outfile, ",\\[Chi]->%.15f", meas_c(M) / h->nv);
  fprintf(outfile, ",\\[Delta]\\[Chi]->%.15f", meas_dc(M) / h->nv);
  fprintf(outfile, ",w->%.15f,\\[Delta]w->%.15f,wc->%.15f,\\[Delta]wc->%.15f,Subscript[n,\"clust\"]->%.15f,Subscript[\\[Delta]n,\"clust\"]->%.15f,Subscript[m,\"clust\"]->%.15f,Subscript[\\[Delta]m,\"clust\"]->%.15f,\\[Tau]->%.15f|>\n", dM->x, meas_dx(dM), meas_c(dM), meas_dc(dM), clust->x / h->nv, meas_dx(clust) / h->nv, meas_c(clust) / h->nv, meas_dc(clust) / h->nv,tau);

  fclose(outfile);

  FILE *image = fopen("out.dat", "a");
  for (v_t i = 0; i < h->nv; i++) {
    fprintf(image, "%" PRIh " ", s->spins[i]);
  }
  fprintf(image, "\n");
  fclose(image);

  free(E);
  free(clust);
  free(H);
  free(s->R);
  free(s->spins);
  graph_free(s->g);
  free(s);
  graph_free(h);
  gsl_rng_free(r);

  return 0;
}

