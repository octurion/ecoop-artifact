#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <emmintrin.h>

#include <array>
#include <vector>

static const double gam = 1.4;
static const double gm1 = gam - 1.0;
static const double gm1i = 1.0 / gm1;
static const double wtg1[2] = {0.5, 0.5};
static const double xi1[2] =
{
	0.211324865405187,
	0.788675134594813,
};
static const double Ng1[4] =
{
	0.788675134594813,
	0.211324865405187,
	0.211324865405187,
	0.788675134594813,
};
static const double Ng1_xi[4] = {-1, -1, 1, 1};
static const double wtg2[4] = {0.25, 0.25, 0.25, 0.25};
static const double Ng2[16] =
{
	0.622008467928146,
	0.166666666666667,
	0.166666666666667,
	0.044658198738520,
	0.166666666666667,
	0.622008467928146,
	0.044658198738520,
	0.166666666666667,
	0.166666666666667,
	0.044658198738520,
	0.622008467928146,
	0.166666666666667,
	0.044658198738520,
	0.166666666666667,
	0.166666666666667,
	0.622008467928146,
};
static const double Ng2_xi[32] = {
	-0.788675134594813,
	 0.788675134594813,
	-0.211324865405187,
	 0.211324865405187,
	-0.788675134594813,
	 0.788675134594813,
	-0.211324865405187,
	 0.211324865405187,
	-0.211324865405187,
	 0.211324865405187,
	-0.788675134594813,
	 0.788675134594813,
	-0.211324865405187,
	 0.211324865405187,
	-0.788675134594813,
	 0.788675134594813,
	-0.788675134594813,
	-0.211324865405187,
	 0.788675134594813,
	 0.211324865405187,
	-0.211324865405187,
	-0.788675134594813,
	 0.211324865405187,
	 0.788675134594813,
	-0.788675134594813,
	-0.211324865405187,
	 0.788675134594813,
	 0.211324865405187,
	-0.211324865405187,
	-0.788675134594813,
	 0.211324865405187,
	 0.788675134594813,
};
static const double minf = 0.1;
static const double m2 = minf * minf;
static const double freq = 1;
static const double kappa = 1;
static const double nmode = 0;
static const double mfan = 1.0;

static const int num_iter = 20;

struct NodePool {
	std::vector<std::array<double, 2>> p_x;
	std::vector<double> p_phim;
	std::vector<double> p_resm;

	std::vector<double> p_V;
	std::vector<double> p_P;
	std::vector<double> p_U;
	size_t size;

	NodePool(size_t num_nodes)
		: p_x(num_nodes)
		, p_phim(num_nodes)
		, p_resm(num_nodes)
		, p_V(num_nodes)
		, p_P(num_nodes)
		, p_U(num_nodes)
		, size(num_nodes)
	{}
};

struct BacknodePool {
	std::vector<size_t> pbedge;  // Node indices

	size_t size;

	BacknodePool(size_t num_backnodes)
		: pbedge(num_backnodes)
		, size(num_backnodes)
	{}
};

struct CellPool {
	std::vector<std::array<size_t, 4>> pcell; // Node indices
	std::vector<std::array<double, 16>> p_K;

	size_t size;

	CellPool(size_t num_cells)
		: pcell(num_cells)
		, p_K(num_cells)
		, size(num_cells)
	{}
};

inline void res_calc(CellPool& cells, NodePool& nodes) {
	for (size_t n = 0; n < cells.size; n++) {
		std::array<const double*, 4> x = {
			nodes.p_x[cells.pcell[n][0]].data(),
			nodes.p_x[cells.pcell[n][1]].data(),
			nodes.p_x[cells.pcell[n][2]].data(),
			nodes.p_x[cells.pcell[n][3]].data(),
		};
		std::array<const double*, 4> phim = {
			&nodes.p_phim[cells.pcell[n][0]],
			&nodes.p_phim[cells.pcell[n][1]],
			&nodes.p_phim[cells.pcell[n][2]],
			&nodes.p_phim[cells.pcell[n][3]],
		};

		std::array<double,16>& K = cells.p_K[n];

		std::array<double*, 4> res = {
			&nodes.p_resm[cells.pcell[n][0]],
			&nodes.p_resm[cells.pcell[n][1]],
			&nodes.p_resm[cells.pcell[n][2]],
			&nodes.p_resm[cells.pcell[n][3]],
		};

		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				K[j * 4 + k] = 0;
			}
		}
		for (int i = 0; i < 4; i++) { // for each gauss point
			double det_x_xi = 0;
			double N_x[8];

			double a = 0;
			for (int m = 0; m < 4; m++)
				det_x_xi += Ng2_xi[4 * i + 16 + m] * x[m][1];
			for (int m = 0; m < 4; m++)
				N_x[m] = det_x_xi * Ng2_xi[4 * i + m];

			a = 0;
			for (int m = 0; m < 4; m++)
				a += Ng2_xi[4 * i + m] * x[m][0];
			for (int m = 0; m < 4; m++)
				N_x[4 + m] = a * Ng2_xi[4 * i + 16 + m];

			det_x_xi *= a;

			a = 0;
			for (int m = 0; m < 4; m++)
				a += Ng2_xi[4 * i + m] * x[m][1];
			for (int m = 0; m < 4; m++)
				N_x[m] -= a * Ng2_xi[4 * i + 16 + m];

			double b = 0;
			for (int m = 0; m < 4; m++)
				b += Ng2_xi[4 * i + 16 + m] * x[m][0];
			for (int m = 0; m < 4; m++)
				N_x[4 + m] -= b * Ng2_xi[4 * i + m];

			det_x_xi -= a * b;

			for (int j = 0; j < 8; j++)
				N_x[j] /= det_x_xi;

			double wt1 = wtg2[i] * det_x_xi;
			// double wt2 = wtg2[i]*det_x_xi/r;

			double u[2] = {0.0, 0.0};
			for (int j = 0; j < 4; j++) {
				u[0] += N_x[j] * phim[j][0];
				u[1] += N_x[4 + j] * phim[j][0];
			}

			double Dk = 1.0 + 0.5 * gm1 * (m2 - (u[0] * u[0] + u[1] * u[1]));
			double rho = pow(Dk, gm1i);
			double rc2 = rho / Dk;

			for (int j = 0; j < 4; j++) {
				res[j][0] += wt1 * rho * (u[0] * N_x[j] + u[1] * N_x[4 + j]);
			}
			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					K[j * 4 + k] +=
						wt1 * rho * (N_x[j] * N_x[k] + N_x[4 + j] * N_x[4 + k]) -
						wt1 * rc2 * (u[0] * N_x[j] + u[1] * N_x[4 + j]) *
						(u[0] * N_x[k] + u[1] * N_x[4 + k]);
				}
			}
		}
	}
}

inline void dirichlet_resm(const BacknodePool& bnodes, NodePool& nodes) {
	for (size_t i = 0; i < bnodes.size; i++) {
		double *res = &nodes.p_resm[bnodes.pbedge[i]];
		*res = 0.0;
	}
}

inline void dirichletPV(const BacknodePool& bnodes, NodePool& nodes) {
	for (size_t i = 0; i < bnodes.size; i++) {
		double *res = &nodes.p_V[bnodes.pbedge[i]];
		*res = 0.0;
	}
}

inline void init_cg(NodePool& nodes, double* c) {
	for (size_t i = 0; i < nodes.size; i++) {
		const double *r = &nodes.p_resm[i];
		double *u = &nodes.p_U[i];
		double *v = &nodes.p_V[i];
		double *p = &nodes.p_P[i];

		*c += (*r) * (*r);
		*p = *r;
		*u = 0;
		*v = 0;
	}
}

inline void spMV(const CellPool& cells, NodePool& nodes) {
	for (size_t i = 0; i < cells.size; i++) {
		std::array<double*, 4> v = {
			&nodes.p_V[cells.pcell[i][0]],
			&nodes.p_V[cells.pcell[i][1]],
			&nodes.p_V[cells.pcell[i][2]],
			&nodes.p_V[cells.pcell[i][3]],
		};
		std::array<const double*, 4> p = {
			&nodes.p_P[cells.pcell[i][0]],
			&nodes.p_P[cells.pcell[i][1]],
			&nodes.p_P[cells.pcell[i][2]],
			&nodes.p_P[cells.pcell[i][3]],
		};

		const std::array<double, 16>& K = cells.p_K[i];

		v[0][0] += K[0] * p[0][0];
		v[0][0] += K[1] * p[1][0];
		v[1][0] += K[1] * p[0][0];
		v[0][0] += K[2] * p[2][0];
		v[2][0] += K[2] * p[0][0];
		v[0][0] += K[3] * p[3][0];
		v[3][0] += K[3] * p[0][0];
		v[1][0] += K[4 + 1] * p[1][0];
		v[1][0] += K[4 + 2] * p[2][0];
		v[2][0] += K[4 + 2] * p[1][0];
		v[1][0] += K[4 + 3] * p[3][0];
		v[3][0] += K[4 + 3] * p[1][0];
		v[2][0] += K[8 + 2] * p[2][0];
		v[2][0] += K[8 + 3] * p[3][0];
		v[3][0] += K[8 + 3] * p[2][0];
		v[3][0] += K[15] * p[3][0];
	}
}

inline void dotPV(const NodePool& nodes, double* c) {
	__m128d acc = _mm_setzero_pd();
	for (size_t i = 0; i < (nodes.size & -2); i += 2) {
		acc += _mm_loadu_pd(&nodes.p_P[i]) * _mm_loadu_pd(&nodes.p_V[i]);
	}
	*c += _mm_cvtsd_f64(acc);
	*c += _mm_cvtsd_f64(_mm_shuffle_pd(acc, acc, 3));

	for (size_t i = (nodes.size & -2); i < nodes.size; i++) {
		const double *p = &nodes.p_P[i];
		const double *v = &nodes.p_V[i];

		*c += (*p) * (*v);
	}
}

inline void updateUR(NodePool& nodes, const double* alpha) {
	auto alpha_pd = _mm_set1_pd(*alpha);
	for (size_t i = 0; i < (nodes.size & -2); i += 2) {
		double *u = &nodes.p_U[i];
		auto u_val = _mm_loadu_pd(u);

		double *r = &nodes.p_resm[i];
		auto r_val = _mm_loadu_pd(r);

		double *v = &nodes.p_V[i];
		auto v_val = _mm_loadu_pd(v);

		const double *p = &nodes.p_P[i];
		auto p_val = _mm_loadu_pd(p);

		_mm_storeu_pd(u, u_val + (alpha_pd * p_val));
		_mm_storeu_pd(r, r_val - (alpha_pd * v_val));
		_mm_storeu_pd(v, _mm_setzero_pd());
	}

	for (size_t i = (nodes.size & -2); i < nodes.size; i++) {
		double *u = &nodes.p_U[i];
		double *r = &nodes.p_resm[i];
		const double *p = &nodes.p_P[i];
		double *v = &nodes.p_V[i];

		*u += (*alpha) * (*p);
		*r -= (*alpha) * (*v);
		*v = 0.0f;
	}
}

inline void dotR(NodePool& nodes, double* c) {
	auto acc = _mm_setzero_pd();
	for (size_t i = 0; i < (nodes.size & -2); i += 2) {
		const double *r = &nodes.p_resm[i];
		auto r_val = _mm_loadu_pd(r);
		
		acc += r_val * r_val;
	}
	*c += _mm_cvtsd_f64(acc);
	*c += _mm_cvtsd_f64(_mm_shuffle_pd(acc, acc, 3));

	for (size_t i = (nodes.size & -2); i < nodes.size; i++) {
		const double *r = &nodes.p_resm[i];

		*c += (*r) * (*r);
	}
}

inline void updateP(NodePool& nodes, const double *beta) {
	auto beta_pd = _mm_set1_pd(*beta);
	for (size_t i = 0; i < (nodes.size & -2); i += 2) {
		const double *r = &nodes.p_resm[i];
		double *p = &nodes.p_P[i];

		_mm_storeu_pd(p, beta_pd * _mm_loadu_pd(p) + _mm_loadu_pd(r));
	}
	for (size_t i = (nodes.size & -2); i < nodes.size; i++) {
		const double *r = &nodes.p_resm[i];
		double *p = &nodes.p_P[i];

		*p = (*beta) * (*p) + (*r);
	}
}

inline void update(NodePool& nodes, double *rms) {
	auto rms_acc = _mm_setzero_pd();
	for (size_t i = 0; i < (nodes.size & -2); i += 2) {
		double *phim = &nodes.p_phim[i];
		auto phim_val = _mm_loadu_pd(phim);

		const double *u = &nodes.p_U[i];
		auto u_val = _mm_loadu_pd(u);

		double *res = &nodes.p_resm[i];

		_mm_storeu_pd(phim, phim_val - u_val);
		_mm_storeu_pd(res, _mm_setzero_pd());
		rms_acc += u_val * u_val;
	}
	*rms += _mm_cvtsd_f64(rms_acc);
	*rms += _mm_cvtsd_f64(_mm_shuffle_pd(rms_acc, rms_acc, 3));
	for (size_t i = (nodes.size & -2); i < nodes.size; i++) {
		double *phim = &nodes.p_phim[i];
		double *res = &nodes.p_resm[i];
		const double *u = &nodes.p_U[i];

		*phim -= *u;
		*res = 0.0;
		*rms += (*u) * (*u);
	}
}

double timespec_elapsed(const struct timespec* end, const struct timespec* start)
{
	return end->tv_sec - start->tv_sec + 1e-9 * (end->tv_nsec - start->tv_nsec);
}

int main(int argc, char** argv)
{
	if (argc != 3) {
		fprintf(stderr, "Usage: %s input_file output_file\n", argv[0]);
		return EXIT_FAILURE;
	}

	FILE* in_file = fopen(argv[1], "r");
	if (in_file == nullptr) {
		perror(argv[1]);
		return EXIT_FAILURE;
	}

	fprintf(stderr, "Now parsing the file...\n");

	size_t nnode, ncell, nbnode;
	if (fscanf(in_file, "%zu %zu %zu \n", &nnode, &ncell, &nbnode) != 3) {
		fprintf(stderr, "Could not parse size of lists of objects\n");
		return EXIT_FAILURE;
	}

	NodePool nodes(nnode);
	BacknodePool bnodes(nbnode);
	CellPool cells(ncell);

	for (size_t n = 0; n < nnode; n++) {
		if (fscanf(in_file, "%lf %lf \n", &nodes.p_x[n][0], &nodes.p_x[n][1]) != 2) {
			fprintf(stderr, "Error reading nodes\n");
			return EXIT_FAILURE;
		}
	}

	for (size_t n = 0; n < ncell; n++) {
		if (fscanf(in_file, "%zu %zu %zu %zu \n",
				   &cells.pcell[n][0],
				   &cells.pcell[n][1],
				   &cells.pcell[n][2],
				   &cells.pcell[n][3]) != 4) {
			fprintf(stderr, "Error reading pcells\n");
			return EXIT_FAILURE;
		}
	}

	for (size_t n = 0; n < nbnode; n++) {
		if (fscanf(in_file, "%zu \n", &bnodes.pbedge[n]) != 1) {
			fprintf(stderr, "Error reading bnodes\n");
			return EXIT_FAILURE;
		}
	}
	fclose(in_file);
	for (size_t i = 0; i < nodes.size; i++) {
		nodes.p_phim[i] = minf * nodes.p_x[i][0];
	}

	fprintf(stderr, "Parsing complete!\n");

	fprintf(stderr, "Now performing iterations...\n");

	struct timespec wall_start, wall_end, cpu_start, cpu_end;

	double res_calc_time = 0;
	double dirichlet_resm_time = 0;
	double init_cg_time = 0;
	double inner_time = 0;
	double spMV_time = 0;

	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpu_start);
	clock_gettime(CLOCK_REALTIME, &wall_start);

	double rms = 1;
	for (int iter = 1; iter <= num_iter; iter++) {
		struct timespec res_calc_start, res_calc_end;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &res_calc_start);
		res_calc(cells, nodes);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &res_calc_end);
		res_calc_time += timespec_elapsed(&res_calc_end, &res_calc_start);

		struct timespec dirichlet_resm_start, dirichlet_resm_end;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &dirichlet_resm_start);
		dirichlet_resm(bnodes, nodes);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &dirichlet_resm_end);
		dirichlet_resm_time += timespec_elapsed(&dirichlet_resm_end, &dirichlet_resm_start);

		double c1 = 0;
		double c2 = 0;
		double c3 = 0;
		double beta = 0;

		struct timespec init_cg_start, init_cg_end;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &init_cg_start);
		init_cg(nodes, &c1);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &init_cg_end);
		init_cg_time += timespec_elapsed(&init_cg_end, &init_cg_start);

		struct timespec inner_start, inner_end;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &inner_start);

		double res0 = sqrt(c1);
		double res = res0;
		int inner_iter = 0;
		int maxiter = 200;
		while (res > 0.1 * res0 && inner_iter < maxiter) {
			struct timespec spMV_start, spMV_end;
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &spMV_start);
			spMV(cells, nodes);
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &spMV_end);
			spMV_time += timespec_elapsed(&spMV_end, &spMV_start);

			dirichletPV(bnodes, nodes);

			c2 = 0;
			dotPV(nodes, &c2);

			double alpha = c1 / c2;
			updateUR(nodes, &alpha);

			c3 = 0;
			dotR(nodes, &c3);

			beta = c3 / c1;
			updateP(nodes, &beta);

			c1 = c3;
			res = sqrt(c1);
			inner_iter++;
		}
		rms = 0;
		update(nodes, &rms);

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &inner_end);
		inner_time += timespec_elapsed(&inner_end, &inner_start);

		printf("rms = %10.5e iter: %d\n", sqrt(rms) / sqrt(nnode), inner_iter);
	}

	clock_gettime(CLOCK_REALTIME, &wall_end);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpu_end);

	fprintf(stderr, "Iterations complete!\n");

	printf("res_calc: %.8f, dirichlet_resm: %.8f,"
			" init_cg: %.8f, spMV: %.8f, inner = %.8f\n",
		res_calc_time, dirichlet_resm_time,
		init_cg_time, spMV_time, inner_time - spMV_time);

	double wall_elapsed = timespec_elapsed(&wall_end, &wall_start);
	printf("Wall clock time: %.8f\n", wall_elapsed);

	double cpu_elapsed = timespec_elapsed(&cpu_end, &cpu_start);
	printf("CPU clock time: %.8f\n", cpu_elapsed);

	fprintf(stderr, "Now writing to file %s...\n", argv[2]);

	FILE* out_file = fopen(argv[2], "w");
	fclose(out_file);
	return EXIT_SUCCESS;
}
