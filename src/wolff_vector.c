
#include <getopt.h>

#include <cluster.h>

double identity(double x) {
  return x;
}

double zero(q_t n, double *H, double *x) {
  return 0.0;
}

double dot(q_t n, double *H, double *x) {
  double total = 0;
  for (q_t i = 0; i < n; i++) {
    total += H[i] * x[i];
  }
  return total;
}

double theta(double x, double y) {
  double val = atan(y / x);

  if (x < 0.0 && y > 0.0) {
    return M_PI + val;
  } else if ( x < 0.0 && y < 0.0 ) {
    return - M_PI + val;
  } else {
    return val;
  }
}

double modulated(q_t n, double *H_info, double *x) {
  return H_info[0] * cos(H_info[1] * theta(x[0], x[1]));
}

double cubic(q_t n, double *H_info, double *x) {
  double tmp = 0;

  for (q_t i = 0; i < n; i++) {
    tmp += pow(x[i], 4);
  }

  return - H_info[0] * tmp;
}

double quadratic(q_t n, double *H_info, double *x) {
  double tmp = 0;

  tmp += pow(x[0], 2);

  for (q_t i = 1; i < n; i++) {
    tmp += - 1.0 / (n - 1.0) * pow(x[i], 2);
  }

  return - 0.5 *  H_info[0] * tmp;
}

int main(int argc, char *argv[]) {

  L_t L = 128;
  count_t N = (count_t)1e7;
  count_t min_runs = 10;
  count_t n = 3;
  q_t q = 2;
  D_t D = 2;
  double T = 2.26918531421;
  double *H = (double *)calloc(MAX_Q, sizeof(double));
  double eps = 0;
  bool silent = false;
  bool record_autocorrelation = false;
  vector_field_t H_type = VECTOR;
  count_t ac_skip = 1;
  count_t W = 10;

  int opt;
  q_t H_ind = 0;

  while ((opt = getopt(argc, argv, "N:n:D:L:q:T:H:m:e:saS:W:f:")) != -1) {
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
    case 'q':
      q = atoi(optarg);
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
    case 'f':
      switch (atoi(optarg)) {
        case 0:
          H_type = VECTOR;
          break;
        case 1:
          H_type = MODULATED;
          break;
        case 2:
          H_type = CUBIC;
          break;
        case 3:
          H_type = QUADRATIC;
          break;
      }
      break;
    default:
      exit(EXIT_FAILURE);
    }
  }

  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, rand_seed());

  vector_state_t *s = (vector_state_t *)calloc(1, sizeof(vector_state_t));

  graph_t *h = graph_create_square(D, L);
  s->g = graph_add_ext(h);

  s->n = q;

  s->spins = (double *)calloc(n * h->nv, sizeof(double));

  for (v_t i = 0; i < h->nv; i++) {
    s->spins[q * i] = 1.0;
  }

  s->H_info = H;
  s->T = T;
  switch (H_type) {
    case VECTOR:
      s->H = dot;
      break;
    case MODULATED:
      s->H = modulated;
      break;
    case CUBIC:
      s->H = cubic;
      break;
    case QUADRATIC:
      s->H = quadratic;
      break;
  }
  s->J = identity;

  s->R = (double *)calloc(q * q, sizeof(double));

  for (q_t i = 0; i < q; i++) {
    s->R[q * i + i] = 1.0;
  }

  s->M = (double *)calloc(q, sizeof(double));
  s->M[0] = 1.0 * (double)h->nv;
  s->E = - ((double)h->ne) * s->J(1.0) - s->H(s->n, s->H_info, s->M);

  double diff = 1e31;
  count_t n_runs = 0;
  count_t n_steps = 0;

  meas_t *E, *clust, **M, *aM;

  M = (meas_t **)malloc(q * sizeof(meas_t *));
  aM = (meas_t *)calloc(q, sizeof(meas_t ));
  for (q_t i = 0; i < q; i++) {
    M[i] = (meas_t *)calloc(1, sizeof(meas_t));
  }

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
           n_runs, fabs(meas_dx(E) / E->x), meas_dx(aM) / aM->x, meas_dc(E) / meas_c(E), meas_dc(aM) / meas_c(aM), h->nv / clust->x);

    count_t n_flips = 0;

    while (n_flips / h->nv < n) {
      v_t v0 = gsl_rng_uniform_int(r, h->nv);
      double *step = gen_rot(r, q);

      v_t tmp_flips = flip_cluster_vector(s, v0, step, r);
      free(step);
      n_flips += tmp_flips;

      if (n_runs > 0) {
        n_steps++;
        meas_update(clust, tmp_flips);

        if (record_autocorrelation && n_steps % ac_skip == 0) {
          update_autocorr(autocorr, s->E);
        }
      }
    }

    double aM_val = 0;

    for (q_t i = 0; i < q; i++) {
      meas_update(M[i], s->M[i]);
      aM_val += s->M[i] * s->M[i];
    }

    meas_update(aM, sqrt(aM_val));
    meas_update(E, s->E);

    diff = fabs(meas_dx(clust) / clust->x);

    n_runs++;
  }

  if (!silent) {
    printf("\033[F\033[J");
  }
  printf("WOLFF: sweep %" PRIu64
         ", dH/H = %.4f, dM/M = %.4f, dC/C = %.4f, dX/X = %.4f, cps: %.1f\n",
         n_runs, fabs(meas_dx(E) / E->x), meas_dx(M[0]) / M[0]->x, meas_dc(E) / meas_c(E), meas_dc(M[0]) / meas_c(M[0]), h->nv / clust->x);

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

    FILE *autocorr_file = fopen("autocorr.dat", "a");

    printf("%g %g\n", Gammas[0], conv_Gamma[0]);

    for (count_t i = 0; i < n+1; i++) {
      fprintf(autocorr_file, "%g ", conv_Gamma[i]);
    }
    fprintf(autocorr_file, "\n");

    fclose(autocorr_file);
    
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

  fprintf(outfile, "<|N->%" PRIcount ",n->%" PRIcount ",D->%" PRID ",L->%" PRIL ",q->%" PRIq ",T->%.15f,H->{", N, n, D, L, q, T);
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "%.15f", H[i]);
    if (i != q-1) {
      fprintf(outfile, ",");
    }
  }
  fprintf(outfile, "},E->%.15f,\\[Delta]E->%.15f,C->%.15f,\\[Delta]C->%.15f,M->{", E->x / h->nv, meas_dx(E) / h->nv, meas_c(E) / h->nv, meas_dc(E) / h->nv);
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "%.15f", M[i]->x / h->nv);
    if (i != q-1) {
      fprintf(outfile, ",");
    }
  }
  fprintf(outfile, "},\\[Delta]M->{");
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "%.15f", meas_dx(M[i]) / h->nv);
    if (i != q-1) {
      fprintf(outfile, ",");
    }
  }
  fprintf(outfile, "},\\[Chi]->{");
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "%.15f", meas_c(M[i]) / h->nv);
    if (i != q-1) {
      fprintf(outfile, ",");
    }
  }
  fprintf(outfile, "},\\[Delta]\\[Chi]->{");
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "%.15f", meas_dc(M[i]) / h->nv);
    if (i != q-1) {
      fprintf(outfile, ",");
    }
  }
  fprintf(outfile, "},aM->%.15f,\\[Delta]aM->%.15f,a\\[Chi]->%.15f,\\[Delta]a\\[Chi]->%.15f,Subscript[n,\"clust\"]->%.15f,Subscript[\\[Delta]n,\"clust\"]->%.15f,Subscript[m,\"clust\"]->%.15f,Subscript[\\[Delta]m,\"clust\"]->%.15f,\\[Tau]->%.15f|>\n", aM->x / h->nv, meas_dx(aM) / h->nv, meas_c(aM) / h->nv, meas_dc(aM) / h->nv, clust->x / h->nv, meas_dx(clust) / h->nv, meas_c(clust) / h->nv, meas_dc(clust) / h->nv,tau);

  fclose(outfile);

  free(E);
  free(clust);
  for (q_t i = 0; i < q; i++) {
    free(M[i]);
  }
  free(M);
  free(H);
  free(s->M);
  free(s->R);
  free(s->spins);
  graph_free(s->g);
  free(s);
  graph_free(h);
  gsl_rng_free(r);

  return 0;
}

