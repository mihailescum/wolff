
#include "wolff.h"

#define TC 2. / log(1. + sqrt(2.))

int main(int argc, char *argv[]) {
	int opt;
	bool use_scho, use_rand_ini;
	lattice_t lat;
	uint16_t L;
	uint32_t n_temps;
	uint64_t N;
	double Tin, Hin, eps, dT, dH;

	L = 128;
	N = 1000;
	lat = SQUARE_LATTICE;
	Tin = 2.3;
	Hin = 0;
	eps = 1e30;
	use_scho = false;
	use_rand_ini = false;
	dT = 0;
	dH = 1;
	n_temps = 1;

	while ((opt = getopt(argc, argv, "N:L:q:T:H:e:St:h:n:r")) != -1) {
		switch (opt) {
			case 'S':
				use_scho = true;
				break;
			case 'r':
				use_rand_ini = true;
				break;
			case 'N':
				N = (uint64_t)atof(optarg);
				break;
			case 'L':
				L = atoi(optarg);
				break;
			case 'T':
				Tin = atof(optarg);
				break;
			case 'H':
				Hin = atof(optarg);
				break;
			case 't':
				dT = atof(optarg);
				break;
			case 'h':
				dH = atof(optarg);
				break;
			case 'n':
				n_temps = atoi(optarg);
				break;
			case 'e':
				eps= atof(optarg);
				break;
			default:
				exit(EXIT_FAILURE);
		}
	}

	gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r, jst_rand_seed());

	graph_t *g = graph_create(lat, TORUS_BOUND, L, false);

	double *Ts = (double *)malloc(n_temps * sizeof(double));
	double *Hs = (double *)malloc(n_temps * sizeof(double));
	double *Tins = (double *)malloc(n_temps * sizeof(double));
	double *Hins = (double *)malloc(n_temps * sizeof(double));

	for (uint32_t i = 0; i < n_temps; i++) {
		double T, H;

		// if Schofield coordinates are used, T and H are treated as R and theta.
		if (use_scho) {
			double R = Tin; double th = Hin;
			double t = R * (1 - pow(th, 2));
			T = TC / (1 - t);
			double h0 = 0.940647;
			double h = h0 * pow(R, 15. / 8.) * hh(th);
			H = h * T;
		} else {
			T = Tin;
			H = Hin;
		}

		Tins[i] = Tin;
		Hins[i] = Hin;

		Ts[i] = T;
		Hs[i] = H;

		Tin += dT;
		if (use_scho) {
			Hin = 1.08144 - (1.08144 - Hin) / dH;
		} else {
			Hin = Hin / dH;
		}
	}

	bool **xs = (bool **)malloc(n_temps * sizeof(bool *));
	for (uint32_t i = 0; i < n_temps; i++) {
		xs[i] = (bool *)calloc(g->nv, sizeof(bool));
	}

	if (use_rand_ini) {
		for (uint32_t j = 0; j < n_temps; j++) {
			for (uint32_t i = 0; i < g->nv; i++) {
				xs[j][i] = gsl_rng_uniform_int(r, 2);
			}
		}
	}

	int32_t *Ms = (int32_t *)calloc(n_temps, sizeof(int32_t));
	for (uint32_t j = 0; j < n_temps; j++) {
		if (use_rand_ini) {
			for (uint32_t i = 0; i < g->nv; i++) {
				if (xs[j][i]) {
					Ms[j]++;
				} else {
					Ms[j]--;
				}
			}
		} else {
			Ms[j] = -g->nv;
		}
	}

	double *Es = (double *)malloc(n_temps * sizeof(double));
	for (uint32_t i = 0; i < n_temps; i++) {
		Es[i] = get_hamiltonian(g, Hs[i], xs[i]);
	}

	double *ps = (double *)malloc(n_temps * sizeof(double));
	for (uint32_t i = 0; i < n_temps; i++) {
		ps[i] = 1 - exp(-2 / Ts[i]);
	}

	double *E1s = (double *)calloc(n_temps, sizeof(double));
	double *E2s = (double *)calloc(n_temps, sizeof(double));
	double *E4s = (double *)calloc(n_temps, sizeof(double));

	double *M1s = (double *)calloc(n_temps, sizeof(double));
	double *M2s = (double *)calloc(n_temps, sizeof(double));
	double *M4s = (double *)calloc(n_temps, sizeof(double));

	double *tE1s = (double *)calloc(n_temps, sizeof(double));
	double *tE2s = (double *)calloc(n_temps, sizeof(double));
	double *tE4s = (double *)calloc(n_temps, sizeof(double));

	double *tM1s = (double *)calloc(n_temps, sizeof(double));
	double *tM2s = (double *)calloc(n_temps, sizeof(double));
	double *tM4s = (double *)calloc(n_temps, sizeof(double));

	uint64_t n_runs = 0;

	double diff = 1e30;

	printf("\n");
	while (diff > eps) {
		printf("\033[F\033[JWOLFF: sweep %llu, %g\n", n_runs, diff);

		for (uint32_t j = 0; j < n_temps; j++) {
			tE1s[j] = 0;
			tE2s[j] = 0;
			tE4s[j] = 0;
			tM1s[j] = 0;
			tM2s[j] = 0;
			tM4s[j] = 0;
		}

		for (uint64_t i = 0; i < N; i++) {
#pragma omp parallel for
			for (uint32_t j = 0; j < n_temps; j++) {
				cluster_t *c = get_cluster(g, xs[j], ps[j], r);
				double dHex = 0;
				double s;

				if (xs[j][c->vi->x]) {
					s = 1;
				} else {
					s = -1;
				}

				dHex = s * Hs[j] * c->nv;

				if (gsl_rng_uniform(r) < exp(-2 * dHex / Ts[j])) {
					while (c->vi != NULL) {
						uint32_t v = queue_del(&(c->vi));
						xs[j][v] = !xs[j][v];
					}
					Es[j] += 2 * c->nb + dHex;
					Ms[j] -= 2 * s * c->nv;
				} else {
					while (c->vi != NULL) {
						queue_del(&(c->vi));
					}
				}

				tE1s[j] += Es[j];
				tM1s[j] += abs(Ms[j]);
				tE2s[j] += pow(Es[j], 2);
				tM2s[j] += pow(Ms[j], 2);
				tE4s[j] += pow(Es[j], 4);
				tM4s[j] += pow(Ms[j], 4);

				free(c);
			}

			for (uint32_t j = 0; j < n_temps - 1; j++) {
				if (gsl_rng_uniform(r) < 0.5) {
					if (gsl_rng_uniform(r) < exp((Es[j + 1] - Es[j]) * (1 / Ts[j + 1] - 1 / Ts[j]))) {
						bool *tmp_x = xs[j];
						double tmp_E = Es[j];
						int32_t tmp_M = Ms[j];

						xs[j] = xs[j + 1];
						xs[j + 1] = tmp_x;

						Es[j] = Es[j + 1];
						Es[j + 1] = tmp_E;

						Ms[j] = Ms[j + 1];
						Ms[j + 1] = tmp_M;
					}
				}
			}
		}

		for (uint32_t j = 0; j < n_temps; j++) {
			tE1s[j] /= N;
			tM1s[j] /= N;
			tE2s[j] /= N;
			tM2s[j] /= N;
			tE4s[j] /= N;
			tM4s[j] /= N;

			if (n_runs > 0) {
				E1s[j] = E1s[j] * ((n_runs - 1.) / n_runs) + tE1s[j] * 1. / n_runs;
				M1s[j] = M1s[j] * ((n_runs - 1.) / n_runs) + tM1s[j] * 1. / n_runs;
				E2s[j] = E2s[j] * ((n_runs - 1.) / n_runs) + tE2s[j] * 1. / n_runs;
				M2s[j] = M2s[j] * ((n_runs - 1.) / n_runs) + tM2s[j] * 1. / n_runs;
				E4s[j] = E4s[j] * ((n_runs - 1.) / n_runs) + tE4s[j] * 1. / n_runs;
				M4s[j] = M4s[j] * ((n_runs - 1.) / n_runs) + tM4s[j] * 1. / n_runs;
			}
		}

		if (n_runs > 0) {
			diff = 0;
			for (uint32_t j = 0; j < n_temps; j++) {
				double dd = sqrt(M2s[j] - pow(M1s[j], 2)) / sqrt(N * n_runs);
				double t_diff = dd / M1s[j];
				if (t_diff > diff) diff = t_diff;
			}
		} else {
			diff = 1e30;
		}

		n_runs++;
	}

	FILE *outfile = fopen("out.dat", "a");

	for (uint32_t j = 0; j < n_temps; j++) {
		double C = (E2s[j] - pow(E1s[j], 2)) / pow(Ts[j], 2);
		double X = (M2s[j] - pow(M1s[j], 2)) / Ts[j];

		double dE1 = sqrt(E2s[j] - pow(E1s[j], 2)) / sqrt(N * n_runs);
		double dM1 = sqrt(M2s[j] - pow(M1s[j], 2)) / sqrt(N * n_runs);
		double dE2 = sqrt(E4s[j] - pow(E2s[j], 2)) / sqrt(N * n_runs);
		double dM2 = sqrt(M4s[j] - pow(M2s[j], 2)) / sqrt(N * n_runs);

		double dC = sqrt(pow(dE2, 2) + pow(2 * E1s[j] * dE1, 2)) / pow(Ts[j], 2);
		double dX = sqrt(pow(dM2, 2) + pow(2 * M1s[j] * dM1, 2)) / Ts[j];

		fprintf(outfile, "%u %u %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n", L, use_scho, Tins[j], Hins[j], E1s[j] / g->nv, dE1 / g->nv, M1s[j] / g->nv, dM1 / g->nv, C / g->nv, dC / g->nv, X / g->nv, dX / g->nv);
	}

	fclose(outfile);

	for (uint32_t i = 0; i < n_temps; i++) {
		free(xs[i]);
	}

	free(xs);

	gsl_rng_free(r);
	graph_free(g);

	free(Ts);
	free(Hs);
	free(Es);
	free(Ms);

	free(E1s);
	free(E2s);
	free(E4s);

	free(M1s);
	free(M2s);
	free(M4s);

	free(tE1s);
	free(tE2s);
	free(tE4s);

	free(tM1s);
	free(tM2s);
	free(tM4s);


	return 0;
}

