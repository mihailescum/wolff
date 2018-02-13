
#include <wolff.h>

int main(int argc, char *argv[]) {

  L_t L = 128;
  count_t N = (count_t)1e7;
  count_t min_runs = 10;
  count_t n = 3;
  q_t q = 2;
  D_t D = 2;
  double T = 2.26918531421;
  double *J = (double *)calloc(MAX_Q, sizeof(double));
  J[0] = 1.0;
  double *H = (double *)calloc(MAX_Q, sizeof(double));
  double eps = 0;
  bool pretend_ising = false;
  bool planar_potts = false;
  bool silent = false;

  int opt;
  q_t J_ind = 0;
  q_t H_ind = 0;

  while ((opt = getopt(argc, argv, "N:n:D:L:q:T:J:H:m:e:Ips")) != -1) {
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
    case 'J':
      J[J_ind] = atof(optarg);
      J_ind++;
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
    case 'I':
      pretend_ising = true;
      break;
    case 'p':
      planar_potts = true;
      break;
    case 's':
      silent = true;
      break;
    default:
      exit(EXIT_FAILURE);
    }
  }

  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, rand_seed());

  if (pretend_ising) {
    q = 2;
    H[1] = -H[0];
    J[1] = -J[0];
  }

  if (planar_potts) {
    for (q_t i = 0; i < q; i++) {
      J[i] = cos(2 * M_PI * i / ((double)q));
    }
  }

  ising_state_t *s = (ising_state_t *)calloc(1, sizeof(ising_state_t));

  graph_t *h = graph_create_square(D, L);
  s->g = graph_add_ext(h);

  s->q = q;

  s->spins = (q_t *)calloc(h->nv + 1, sizeof(q_t));

  s->T = T;
  s->H = H;
  s->J = J;

  s->J_probs = (double *)calloc(pow(q, 2), sizeof(double));
  for (q_t i = 0; i < q; i++) {
    for (q_t j = 0; j < q; j++) {
      s->J_probs[q * i + j] = 1.0 - exp((s->J[i] - s->J[j]) / T);
    }
  }
  s->H_probs = (double *)calloc(pow(q, 2), sizeof(double));
  for (q_t i = 0; i < q; i++) {
    for (q_t j = 0; j < q; j++) {
      s->H_probs[q * i + j] = 1.0 - exp((s->H[i] - s->H[j]) / T);
    }
  }

  s->M = (v_t *)calloc(q, sizeof(v_t));
  s->M[0] = h->nv;
  s->E = - ((double)h->ne) * s->J[0] - ((double)h->nv) * s->H[0];

  double diff = 1e31;
  count_t n_runs = 0;

  meas_t *E, *clust, **M, **sE, ***sM, **lifetimes;

  M = (meas_t **)malloc(q * sizeof(meas_t *));
  lifetimes = (meas_t **)malloc(q * sizeof(meas_t *));
  for (q_t i = 0; i < q; i++) {
    M[i] = (meas_t *)calloc(1, sizeof(meas_t));
    lifetimes[i] = (meas_t *)calloc(1, sizeof(meas_t));
  }

  E = calloc(1, sizeof(meas_t));
  clust = calloc(1, sizeof(meas_t));

  sE = (meas_t **)malloc(q * sizeof(meas_t *));
  sM = (meas_t ***)malloc(q * sizeof(meas_t **));

  for (q_t i = 0; i < q; i++) {
    sE[i] = (meas_t *)calloc(1, sizeof(meas_t));
    sM[i] = (meas_t **)malloc(q * sizeof(meas_t *));
    for (q_t j = 0; j < q; j++) {
      sM[i][j] = (meas_t *)calloc(1, sizeof(meas_t));
    }
  }

  count_t *freqs = (count_t *)calloc(q, sizeof(count_t));
  count_t lifetime_n = 0;
  q_t cur_M = 0;

  if (!silent) printf("\n");
  while (((diff > eps || diff != diff) && n_runs < N) || n_runs < min_runs) {
    if (!silent) printf("\033[F\033[JWOLFF: sweep %" PRIu64
           ", dH/H = %.4f, dM/M = %.4f, dC/C = %.4f, dX/X = %.4f, cps: %.1f\n",
           n_runs, fabs(E->dx / E->x), M[0]->dx / M[0]->x, E->dc / E->c, M[0]->dc / M[0]->c, h->nv / clust->x);

    count_t n_flips = 0;

    q_t max_M_i;
    v_t max_M;
    q_t n_at_max;

    while (n_flips / h->nv < n) {
      v_t v0 = gsl_rng_uniform_int(r, h->nv);
      q_t step = 1 + gsl_rng_uniform_int(r, q - 1);

      v_t tmp_flips = flip_cluster(s, v0, step, r);
      n_flips += tmp_flips;

      max_M_i = 0;
      max_M = 0;
      n_at_max = 0;

      for (q_t i = 0; i < q; i++) {
        if (s->M[i] > max_M) {
          max_M = s->M[i];
          max_M_i = i;
          n_at_max = 1;
        } else if (s->M[i] == max_M) {
          n_at_max++;
        }
      }

      if (n_at_max == 1) {
        if (max_M_i == cur_M) {
          lifetime_n++;
        } else {
          if (cur_M != MAX_Q) {
            update_meas(lifetimes[cur_M], lifetime_n);
          }
          lifetime_n = 0;
          cur_M = max_M_i;
        }
      } else {
        if (cur_M != MAX_Q) {
          update_meas(lifetimes[cur_M], lifetime_n);
          cur_M = MAX_Q;
        }
      }

      update_meas(clust, tmp_flips);
    }

    for (q_t i = 0; i < q; i++) {
      update_meas(M[i], s->M[i]);
    }
    update_meas(E, s->E);

    if (n_at_max == 1) {
      for (q_t i = 0; i < q; i++) {
        update_meas(sM[max_M_i][i], s->M[i]);
      }
      update_meas(sE[max_M_i], s->E);
      freqs[max_M_i]++;
    }

    diff = fabs(sM[0][0]->dc / sM[0][0]->c);

    n_runs++;
  }
  if (!silent) {
    printf("\033[F\033[J");
  }
  printf("WOLFF: sweep %" PRIu64
         ", dH/H = %.4f, dM/M = %.4f, dC/C = %.4f, dX/X = %.4f, cps: %.1f\n",
         n_runs, fabs(E->dx / E->x), M[0]->dx / M[0]->x, E->dc / E->c, M[0]->dc / M[0]->c, h->nv / clust->x);

  FILE *outfile = fopen("out.m", "a");

  fprintf(outfile, "<|D->%" PRID ",L->%" PRIL ",q->%" PRIq ",T->%.15f,J->{", D, L, q, T);
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "%.15f", J[i]);
    if (i != q-1) {
      fprintf(outfile, ",");
    }
  }
  fprintf(outfile, "},H->{");
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "%.15f", H[i]);
    if (i != q-1) {
      fprintf(outfile, ",");
    }
  }
  fprintf(outfile, "},E->%.15f,\\[Delta]E->%.15f,C->%.15f,\\[Delta]C->%.15f,M->{", E->x / h->nv, E->dx / h->nv, E->c / h->nv, E->dc / h->nv);
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "%.15f", M[i]->x / h->nv);
    if (i != q-1) {
      fprintf(outfile, ",");
    }
  }
  fprintf(outfile, "},\\[Delta]M->{");
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "%.15f", M[i]->dx / h->nv);
    if (i != q-1) {
      fprintf(outfile, ",");
    }
  }
  fprintf(outfile, "},\\[Chi]->{");
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "%.15f", M[i]->c / h->nv);
    if (i != q-1) {
      fprintf(outfile, ",");
    }
  }
  fprintf(outfile, "},\\[Delta]\\[Chi]->{");
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "%.15f", M[i]->dc / h->nv);
    if (i != q-1) {
      fprintf(outfile, ",");
    }
  }
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "},Subscript[E,%" PRIq "]->%.15f,Subscript[\\[Delta]E,%" PRIq "]->%.15f,Subscript[C,%" PRIq "]->%.15f,Subscript[\\[Delta]C,%" PRIq "]->%.15f,Subscript[M,%" PRIq "]->{", i, sE[i]->x / h->nv, i, sE[i]->dx / h->nv, i, sE[i]->c / h->nv, i, sE[i]->dc / h->nv, i);
    for (q_t j = 0; j < q; j++) {
      fprintf(outfile, "%.15f", sM[i][j]->x / h->nv);
      if (j != q-1) {
        fprintf(outfile, ",");
      }
    }
    fprintf(outfile, "},Subscript[\\[Delta]M,%" PRIq "]->{", i);
    for (q_t j = 0; j < q; j++) {
      fprintf(outfile, "%.15f", sM[i][j]->dx / h->nv);
      if (j != q-1) {
        fprintf(outfile, ",");
      }
    }
    fprintf(outfile, "},Subscript[\\[Chi],%" PRIq "]->{", i);
    for (q_t j = 0; j < q; j++) {
      fprintf(outfile, "%.15f", sM[i][j]->c / h->nv);
      if (j != q-1) {
        fprintf(outfile, ",");
      }
    }
    fprintf(outfile, "},Subscript[\\[Delta]\\[Chi],%" PRIq "]->{", i);
    for (q_t j = 0; j < q; j++) {
      fprintf(outfile, "%.15f", sM[i][j]->dc / h->nv);
      if (j != q-1) {
        fprintf(outfile, ",");
      }
    }
  }
  fprintf(outfile,"}");
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, ",Subscript[f,%" PRIq "]->%.15f,Subscript[\\[Delta]f,%" PRIq "]->%.15f", i, (double)freqs[i] / (double)n_runs, i, sqrt(freqs[i]) / (double)n_runs);
  }
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, ",Subscript[t,%" PRIq "]->%.15f,Subscript[\\[Delta]t,%" PRIq "]->%.15f", i, lifetimes[i]->x, i, lifetimes[i]->dx);
  }
  fprintf(outfile, ",Subscript[n,\"clust\"]->%.15f,Subscript[\\[Delta]n,\"clust\"]->%.15f|>\n", clust->x / h->nv, clust->dx / h->nv);

  fclose(outfile);

  free(E);
  free(clust);
  for (q_t i = 0; i < q; i++) {
    free(M[i]);
    for (q_t j = 0; j < q; j++) {
      free(sM[i][j]);
    }
    free(sM[i]);
  }
  free(M);
  free(sM);
  for (q_t i = 0; i < q; i++) {
    free(sE[i]);
    free(lifetimes[i]);
  }
  free(freqs);
  free(sE);
  free(s->H_probs);
  free(s->J_probs);
  free(s->M);
  free(s->spins);
  graph_free(s->g);
  free(s);
  free(H);
  free(J);
  graph_free(h);
  gsl_rng_free(r);

  return 0;
}

