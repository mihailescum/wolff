
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

  int opt;
  q_t J_ind = 0;
  q_t H_ind = 0;

  while ((opt = getopt(argc, argv, "N:n:D:L:q:T:J:H:m:e:Ip")) != -1) {
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

  meas_t *E, *clust, **M, **sE, ***sM;

  M = (meas_t **)malloc(q * sizeof(meas_t *));
  for (q_t i = 0; i < q; i++) {
    M[i] = (meas_t *)calloc(1, sizeof(meas_t));
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

  printf("\n");
  while (((diff > eps || diff != diff) && n_runs < N) || n_runs < min_runs) {
    printf("\033[F\033[JWOLFF: sweep %" PRIu64
           ", dH/H = %.4f, dM/M = %.4f, dC/C = %.4f, dX/X = %.4f, cps: %.1f\n",
           n_runs, fabs(E->dx / E->x), M[0]->dx / M[0]->x, E->dc / E->c, M[0]->dc / M[0]->c, h->nv / clust->x);

    count_t n_flips = 0;

    while (n_flips / h->nv < n) {
      v_t v0 = gsl_rng_uniform_int(r, h->nv);
      q_t step = 1 + gsl_rng_uniform_int(r, q - 1);

      v_t tmp_flips = flip_cluster(s, v0, step, r);
      n_flips += tmp_flips;

      update_meas(clust, tmp_flips);
    }

    for (q_t i = 0; i < q; i++) {
      update_meas(M[i], s->M[i]);
    }
    update_meas(E, s->E);

    q_t max_M_i = 0;
    double max_M = 0;

    for (q_t i = 0; i < q; i++) {
      if (s->M[i] > max_M) {
        max_M = s->M[i];
        max_M_i = i;
      }
    }

    for (q_t i = 0; i < q; i++) {
      update_meas(sM[max_M_i][i], s->M[i]);
    }
    update_meas(sE[max_M_i], s->E);

    diff = fabs(sM[0][0]->dc / sM[0][0]->c);

    n_runs++;
  }

  printf("\033[F\033[JWOLFF: sweep %" PRIu64
         ", dH/H = %.4f, dM/M = %.4f, dC/C = %.4f, dX/X = %.4f, cps: %.1f\n",
         n_runs, fabs(E->dx / E->x), M[0]->dx / M[0]->x, E->dc / E->c, M[0]->dc / M[0]->c, h->nv / clust->x);

  FILE *outfile = fopen("out.m", "a");

  fprintf(outfile, "{\"D\"->%" PRID ",\"L\"->%" PRIL ",\"q\"->%" PRIq ",\"T\"->%.15f,\"J\"->{", D, L, q, T);
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "%.15f", J[i]);
    if (i != q-1) {
      fprintf(outfile, ",");
    }
  }
  fprintf(outfile, "},\"H\"->{");
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "%.15f", H[i]);
    if (i != q-1) {
      fprintf(outfile, ",");
    }
  }
  fprintf(outfile, "},\"E\"->{%.15f,%.15f},\"C\"->{%.15f,%.15f},\"M\"->{", E->x / h->nv, E->dx / h->nv, E->c / h->nv, E->dc / h->nv);
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "{%.15f,%.15f}", M[i]->x / h->nv, M[i]->dx / h->nv);
    if (i != q-1) {
      fprintf(outfile, ",");
    }
  }
  fprintf(outfile, "},\"\\[Chi]\"->{");
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "{%.15f,%.15f}", M[i]->c / h->nv, M[i]->dc / h->nv);
    if (i != q-1) {
      fprintf(outfile, ",");
    }
  }
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "},\"sE%" PRIq "\"->{%.15f,%.15f},\"C%" PRIq "\"->{%.15f,%.15f},\"M%" PRIq "\"->{", i, sE[i]->x / h->nv, sE[i]->dx / h->nv, i, sE[i]->c / h->nv, sE[i]->dc / h->nv, i);
    for (q_t j = 0; j < q; j++) {
      fprintf(outfile, "{%.15f,%.15f}", sM[i][j]->x / h->nv, sM[i][j]->dx / h->nv);
      if (j != q-1) {
        fprintf(outfile, ",");
      }
    }
    fprintf(outfile, "},\"\\[Chi]%" PRIq "\"->{", i);
    for (q_t j = 0; j < q; j++) {
      fprintf(outfile, "{%.15f,%.15f}", sM[i][j]->c / h->nv, sM[i][j]->dc / h->nv);
      if (j != q-1) {
        fprintf(outfile, ",");
      }
    }
  }
  fprintf(outfile, "},\"n\"->{%.15f,%.15f}}\n", clust->c / h->nv, clust->dc / h->nv);

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
  }
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

