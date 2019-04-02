#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <algorithm>
#include <array>
#include <unordered_set>
#include <vector>

static const size_t CHUNK_SIZE = 384;

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

struct CellRaw {
	std::array<size_t, 4> pcell;
};

struct Node {
	std::array<double, 2> p_x;
	double p_phim = 0;
	double p_resm = 0;
	double p_V = 0;
	double p_P = 0;
	double p_U = 0;
};

struct Backnode {
	size_t pbedge = 0;  // Node indices
};

struct Cell {
	std::array<size_t, 4> pcell;
	std::array<double, 16> p_K;
};

struct CellSet {
	std::vector<size_t> indices;
	std::unordered_set<size_t> set;
};

struct Color {
	std::vector<std::pair<size_t, size_t>> chunk_ranges;
	std::unordered_set<size_t> indices;
};

struct Chunk {
	std::pair<size_t, size_t> range;
	std::vector<size_t> indices;
};

void res_calc(std::vector<Cell>& cells, std::vector<Node>& nodes, const std::vector<std::pair<size_t, size_t>>& ranges) {
	#pragma omp parallel
	for (const auto& r: ranges) {
		size_t begin = r.first;
		size_t end = r.second;

		#pragma omp for
		for (size_t n = begin; n < end; n++) {
			std::array<const double*, 4> x = {
				nodes[cells[n].pcell[0]].p_x.data(),
				nodes[cells[n].pcell[1]].p_x.data(),
				nodes[cells[n].pcell[2]].p_x.data(),
				nodes[cells[n].pcell[3]].p_x.data(),
			};
			std::array<const double*, 4> phim = {
				&nodes[cells[n].pcell[0]].p_phim,
				&nodes[cells[n].pcell[1]].p_phim,
				&nodes[cells[n].pcell[2]].p_phim,
				&nodes[cells[n].pcell[3]].p_phim,
			};

			std::array<double,16>& K = cells[n].p_K;

			std::array<double*, 4> res = {
				&nodes[cells[n].pcell[0]].p_resm,
				&nodes[cells[n].pcell[1]].p_resm,
				&nodes[cells[n].pcell[2]].p_resm,
				&nodes[cells[n].pcell[3]].p_resm,
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
}

void dirichlet_resm(const std::vector<Backnode>& bnodes, std::vector<Node>& nodes) {
	for (size_t i = 0; i < bnodes.size(); i++) {
		double *res = &nodes[bnodes[i].pbedge].p_resm;
		*res = 0.0;
	}
}

void dirichletPV(const std::vector<Backnode>& bnodes, std::vector<Node>& nodes) {
	for (size_t i = 0; i < bnodes.size(); i++) {
		double *res = &nodes[bnodes[i].pbedge].p_V;
		*res = 0.0;
	}
}

void init_cg(std::vector<Node>& nodes, double* c) {
	double c_acc = 0;

	#pragma omp parallel for reduction(+:c_acc)
	for (size_t i = 0; i < nodes.size(); i++) {
		const double *r = &nodes[i].p_resm;
		double *u = &nodes[i].p_U;
		double *v = &nodes[i].p_V;
		double *p = &nodes[i].p_P;

		c_acc += (*r) * (*r);
		*p = *r;
		*u = 0;
		*v = 0;
	}

	*c += c_acc;
}

void spMV(const std::vector<Cell>& cells, std::vector<Node>& nodes, const std::vector<std::pair<size_t, size_t>>& ranges) {
	#pragma omp parallel
	for (const auto& r: ranges) {
		size_t begin = r.first;
		size_t end = r.second;

		#pragma omp for
		for (size_t i = begin; i < end; i++) {
			std::array<double*, 4> v = {
				&nodes[cells[i].pcell[0]].p_V,
				&nodes[cells[i].pcell[1]].p_V,
				&nodes[cells[i].pcell[2]].p_V,
				&nodes[cells[i].pcell[3]].p_V,
			};
			std::array<const double*, 4> p = {
				&nodes[cells[i].pcell[0]].p_P,
				&nodes[cells[i].pcell[1]].p_P,
				&nodes[cells[i].pcell[2]].p_P,
				&nodes[cells[i].pcell[3]].p_P,
			};

			const std::array<double, 16>& K = cells[i].p_K;

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
}

void dotPV(const std::vector<Node>& nodes, double* c) {
	double c_acc = 0;

	#pragma omp parallel for reduction(+:c_acc)
	for (size_t i = 0; i < nodes.size(); i++) {
		const double *p = &nodes[i].p_P;
		const double *v = &nodes[i].p_V;

		c_acc += (*p) * (*v);
	}

	*c += c_acc;
}

void updateUR(std::vector<Node>& nodes, const double* alpha) {
	#pragma omp parallel for
	for (size_t i = 0; i < nodes.size(); i++) {
		double *u = &nodes[i].p_U;
		double *r = &nodes[i].p_resm;
		const double *p = &nodes[i].p_P;
		double *v = &nodes[i].p_V;

		*u += (*alpha) * (*p);
		*r -= (*alpha) * (*v);
		*v = 0.0f;
	}
}

void dotR(std::vector<Node>& nodes, double* c) {
	double c_acc = 0;

	#pragma omp parallel for reduction(+:c_acc)
	for (size_t i = 0; i < nodes.size(); i++) {
		const double *r = &nodes[i].p_resm;

		c_acc += (*r) * (*r);
	}

	*c += c_acc;
}

void updateP(std::vector<Node>& nodes, const double *beta) {
	#pragma omp parallel for
	for (size_t i = 0; i < nodes.size(); i++) {
		const double *r = &nodes[i].p_resm;
		double *p = &nodes[i].p_P;

		*p = (*beta) * (*p) + (*r);
	}
}

void update(std::vector<Node>& nodes, double *rms) {
	double rms_acc = 0;

	#pragma omp parallel for reduction(+:rms_acc)
	for (size_t i = 0; i < nodes.size(); i++) {
		double *phim = &nodes[i].p_phim;
		double *res = &nodes[i].p_resm;
		const double *u = &nodes[i].p_U;

		*phim -= *u;
		*res = 0.0;
		rms_acc += (*u) * (*u);
	}

	*rms += rms_acc;
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

	std::vector<Node> nodes(nnode);
	std::vector<Backnode> bnodes(nbnode);
	std::vector<Cell> cells(ncell);

	for (size_t n = 0; n < nnode; n++) {
		if (fscanf(in_file, "%lf %lf \n", &nodes[n].p_x[0], &nodes[n].p_x[1]) != 2) {
			fprintf(stderr, "Error reading nodes\n");
			return EXIT_FAILURE;
		}
	}

	std::vector<CellRaw> cells_raw(ncell);
	for (size_t n = 0; n < ncell; n++) {
		if (fscanf(in_file, "%zu %zu %zu %zu \n",
				   &cells_raw[n].pcell[0],
				   &cells_raw[n].pcell[1],
				   &cells_raw[n].pcell[2],
				   &cells_raw[n].pcell[3]) != 4) {
			fprintf(stderr, "Error reading pcells\n");
			return EXIT_FAILURE;
		}
	}

	std::vector<std::pair<size_t, size_t>> ranges;
	std::vector<Color> colors;

	{
		size_t idx = 0;
		for (size_t chunk_id = 0; chunk_id < (ncell + CHUNK_SIZE - 1) / CHUNK_SIZE; chunk_id++) {
			Chunk chunk;
			chunk.range.first = idx;
			chunk.indices.reserve(2 * CHUNK_SIZE);
			for (size_t block_idx = 0; idx < ncell && block_idx < CHUNK_SIZE; idx++, block_idx++) {
				for (size_t k = 0; k < 4; k++) {
					chunk.indices.push_back(cells_raw[idx].pcell[k]);
				}
			}
			chunk.range.second = idx;
			std::sort(chunk.indices.begin(), chunk.indices.end());
			chunk.indices.erase(std::unique(chunk.indices.begin(), chunk.indices.end()), chunk.indices.end());

			Color* found = nullptr;
			for (auto& c: colors) {
				bool intersects = std::any_of(chunk.indices.begin(), chunk.indices.end(), [&c](size_t k) {
					return c.indices.find(k) != c.indices.end();
				});
				if (!intersects) {
					found = &c;
					break;
				}
			}

			if (found == nullptr) {
				colors.emplace_back();
				found = &colors.back();
			}

			found->chunk_ranges.push_back(chunk.range);
			found->indices.insert(chunk.indices.begin(), chunk.indices.end());
		}
	}

	{
		size_t last_idx = 0;
		size_t i = 0;
		for (const auto& c: colors) {
			for (const auto& r: c.chunk_ranges) {
				for (size_t k = r.first; k < r.second; k++, i++) {
					for (size_t l = 0; l < 4; l++) {
						cells[i].pcell[l] = cells_raw[k].pcell[l];
					}
				}
			}

			ranges.emplace_back(last_idx, i);
			last_idx = i;
		}
	}

	for (size_t n = 0; n < nbnode; n++) {
		if (fscanf(in_file, "%zu \n", &bnodes[n].pbedge) != 1) {
			fprintf(stderr, "Error reading bnodes\n");
			return EXIT_FAILURE;
		}
	}
	fclose(in_file);
	for (size_t i = 0; i < nodes.size(); i++) {
		nodes[i].p_phim = minf * nodes[i].p_x[0];
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
		res_calc(cells, nodes, ranges);
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
			spMV(cells, nodes, ranges);
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
