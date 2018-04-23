
#include <getopt.h>

#include <cluster.h>

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
  bool snapshots = false;
  bool snapshot = false;
  bool record_autocorrelation = false;
  bool record_distribution = false;
  count_t W = 10;
  count_t ac_skip = 1;

  int opt;
  q_t J_ind = 0;
  q_t H_ind = 0;

  while ((opt = getopt(argc, argv, "N:n:D:L:q:T:J:H:m:e:IpsSPak:W:d")) != -1) {
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
    case 'S':
      snapshots = true;
      break;
    case 'P':
      snapshot = true;
      break;
    case 'a':
      record_autocorrelation = true;
      break;
    case 'k':
      ac_skip = (count_t)atof(optarg);
      break;
    case 'W':
      W = (count_t)atof(optarg);
      break;
    case 'd':
      record_distribution = true;
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

  s->spins = (q_t *)calloc(h->nv, sizeof(q_t));

  s->T = T;
  s->H = H;
  s->J = J;
  s->R = (dihedral_t *)calloc(1, sizeof(dihedral_t));

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
  count_t n_steps = 0;

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

  count_t *freqs = (count_t *)calloc(q, sizeof(count_t));
  q_t cur_M = 0;

  autocorr_t *autocorr;
  if (record_autocorrelation) {
    autocorr = (autocorr_t *)calloc(1, sizeof(autocorr_t));
    autocorr->W = 2 * W + 1;
    autocorr->OO = (double *)calloc(2 * W + 1, sizeof(double));
  }

  count_t *cluster_dist;
  if (record_distribution) {
    cluster_dist = (count_t *)calloc(h->nv, sizeof(count_t));
  }

  if (!silent) printf("\n");
  while (((diff > eps || diff != diff) && n_runs < N) || n_runs < min_runs) {
    if (!silent) printf("\033[F\033[JWOLFF: sweep %" PRIu64
           ", dH/H = %.4f, dM/M = %.4f, dC/C = %.4f, dX/X = %.4f, cps: %.1f\n",
           n_runs, fabs(meas_dx(E) / E->x), meas_dx(M[0]) / M[0]->x, meas_dc(E) / meas_c(E), meas_dc(M[0]) / meas_c(M[0]), h->nv / clust->x);

    count_t n_flips = 0;

    while (n_flips / h->nv < n) {
      v_t v0 = gsl_rng_uniform_int(r, h->nv);
      q_t step;
     
      if (q == 2) {
        step = 1;
      } else {
        step = gsl_rng_uniform_int(r, q);
      }

      v_t tmp_flips = flip_cluster(s, v0, step, r);
      n_flips += tmp_flips;

      if (n_runs > 0) {
        n_steps++;
        meas_update(clust, tmp_flips);

        if (record_autocorrelation && n_steps % ac_skip == 0) {
          update_autocorr(autocorr, s->E);
        }

        if (record_distribution) {
          cluster_dist[tmp_flips - 1]++;
        }
      }

    }

    for (q_t i = 0; i < q; i++) {
      meas_update(M[i], s->M[i]);
    }
    meas_update(E, s->E);

    q_t n_at_max = 0;
    q_t max_M_i = 0;
    v_t max_M = 0;

    for (q_t i = 0; i < q; i++) {
      if (s->M[i] > max_M) {
        n_at_max = 1;
        max_M_i = i;
        max_M = s->M[i];
      } else if (s->M[i] == max_M) {
        n_at_max++;
      }
    }

    if (n_at_max == 1) {
      for (q_t i = 0; i < q; i++) {
        meas_update(sM[max_M_i][i], s->M[i]);
      }
      meas_update(sE[max_M_i], s->E);
      freqs[max_M_i]++;
    }

    diff = fabs(meas_dc(sM[0][0]) / meas_c(sM[0][0]));

    n_runs++;
  }
  if (!silent) {
    printf("\033[F\033[J");
  }
  printf("WOLFF: sweep %" PRIu64
         ", dH/H = %.4f, dM/M = %.4f, dC/C = %.4f, dX/X = %.4f, cps: %.1f\n",
         n_runs, fabs(meas_dx(E) / E->x), meas_dx(M[0]) / M[0]->x, meas_dc(E) / meas_c(E), meas_dc(M[0]) / meas_c(M[0]), h->nv / clust->x);

  if (snapshots) {
    FILE *snapfile = fopen("snapshots.m", "a");
    fprintf(snapfile, "\n");
  }

  if (snapshot) {
    FILE *snapfile = fopen("snapshot.m", "a");
    fprintf(snapfile, "{{");
    for (L_t i = 0; i < L; i++) {
      fprintf(snapfile, "{");
      for (L_t j = 0; j < L; j++) {
        fprintf(snapfile, "%" PRIq, dihedral_inverse_act(q, s->R, s->spins[L * i + j]));
        if (j != L - 1) {
          fprintf(snapfile, ",");
        }
      }
      fprintf(snapfile, "}");
      if (i != L - 1) {
        fprintf(snapfile, ",");
      }
    }
    fprintf(snapfile, "},{%" PRIq ",%d}}\n", s->R->i, s->R->r);
    fclose(snapfile);
  }

  double tau = 0;
  double tauyw;
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
    
    tau = ttau * ac_skip * clust->x / h->nv;
    tauyw = yule_walker(autocorr) * ac_skip * clust->x / h->nv;
    
    free(Gammas);
    free(autocorr->OO);
    while (autocorr->Op != NULL) {
      stack_pop_d(&(autocorr->Op));
    }
    free(autocorr);
  }

  if (tau_failed) {
    tau = 0;
  }

  FILE *outfile = fopen("out.m", "a");

  fprintf(outfile, "<|N->%" PRIcount ",n->%" PRIcount ",D->%" PRID ",L->%" PRIL ",q->%" PRIq ",T->%.15f,J->{", N, n, D, L, q, T);
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
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "},Subscript[E,%" PRIq "]->%.15f,Subscript[\\[Delta]E,%" PRIq "]->%.15f,Subscript[C,%" PRIq "]->%.15f,Subscript[\\[Delta]C,%" PRIq "]->%.15f,Subscript[M,%" PRIq "]->{", i, sE[i]->x / h->nv, i, meas_dx(sE[i]) / h->nv, i, meas_c(sE[i]) / h->nv, i, meas_dc(sE[i]) / h->nv, i);
    for (q_t j = 0; j < q; j++) {
      fprintf(outfile, "%.15f", sM[i][j]->x / h->nv);
      if (j != q-1) {
        fprintf(outfile, ",");
      }
    }
    fprintf(outfile, "},Subscript[\\[Delta]M,%" PRIq "]->{", i);
    for (q_t j = 0; j < q; j++) {
      fprintf(outfile, "%.15f", meas_dx(sM[i][j]) / h->nv);
      if (j != q-1) {
        fprintf(outfile, ",");
      }
    }
    fprintf(outfile, "},Subscript[\\[Chi],%" PRIq "]->{", i);
    for (q_t j = 0; j < q; j++) {
      fprintf(outfile, "%.15f", meas_c(sM[i][j]) / h->nv);
      if (j != q-1) {
        fprintf(outfile, ",");
      }
    }
    fprintf(outfile, "},Subscript[\\[Delta]\\[Chi],%" PRIq "]->{", i);
    for (q_t j = 0; j < q; j++) {
      fprintf(outfile, "%.15f", meas_dc(sM[i][j]) / h->nv);
      if (j != q-1) {
        fprintf(outfile, ",");
      }
    }
  }
  fprintf(outfile,"}");
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, ",Subscript[f,%" PRIq "]->%.15f,Subscript[\\[Delta]f,%" PRIq "]->%.15f", i, (double)freqs[i] / (double)n_runs, i, sqrt(freqs[i]) / (double)n_runs);
  }
  fprintf(outfile, ",Subscript[n,\"clust\"]->%.15f,Subscript[\\[Delta]n,\"clust\"]->%.15f,Subscript[m,\"clust\"]->%.15f,Subscript[\\[Delta]m,\"clust\"]->%.15f,\\[Tau]->%.15f,\\[Tau]yw->%.15f", clust->x / h->nv, meas_dx(clust) / h->nv, meas_c(clust) / h->nv, meas_dc(clust) / h->nv,tau, tauyw);
  if (record_distribution) {
    fprintf(outfile, ",S->{");
    for (v_t i = 0; i < h->nv; i++) {
      fprintf(outfile, "%" PRIcount, cluster_dist[i]);
      if (i != h->nv - 1) {
        fprintf(outfile, ",");
      }
    }
    fprintf(outfile, "}");
    free(cluster_dist);
  }
  fprintf(outfile, "|>\n");

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
  free(freqs);
  free(sE);
  free(s->H_probs);
  free(s->J_probs);
  free(s->M);
  free(s->spins);
  free(s->R);
  graph_free(s->g);
  free(s);
  free(H);
  free(J);
  graph_free(h);
  gsl_rng_free(r);

  return 0;
}

