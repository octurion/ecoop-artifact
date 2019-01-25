#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <array>
#include <vector>

static const int num_iter = 1000;

static const double gam = 1.4f;
static const double gm1 = gam - 1.0f;
static const double cfl = 0.9f;
static const double eps = 0.05f;

static const double mach = 0.4f;
static const double alpha = 3.0f * atan(1.0f) / 45.0f;
static const double p = 1.0f;
static const double r = 1.0f;
static const double u = sqrt(gam * p / r) * mach;
static const double e = p / (r * gm1) + 0.5f * u * u;

static const double qinf[4] = {
	r, r * u, 0.0f, r * e,
};

struct EdgePool {
	std::vector<std::array<size_t, 2>> pedge;  // Node indices
	std::vector<std::array<size_t, 2>> pecell; // Cell indices
	size_t size;

	EdgePool(size_t num_edges)
		: pedge(num_edges)
		, pecell(num_edges)
		, size(num_edges)
	{}
};

struct BackedgePool {
	std::vector<std::array<size_t, 2>> pbedge; // Node indices
	std::vector<size_t> pbecell;               // Cell indices

	std::vector<int> p_bound; // Bound
	size_t size;

	BackedgePool(size_t num_backedges)
		: pbedge(num_backedges)
		, pbecell(num_backedges)
		, p_bound(num_backedges)
		, size(num_backedges)
	{}
};

struct CellPool {
	std::vector<std::array<size_t, 4>> pcell; // Node indices
	std::vector<std::array<std::array<double, 2>, 4>> pcell_new; // Node indices

	std::vector<std::array<double, 4>> p_q;
	std::vector<std::array<double, 4>> qold;
	std::vector<double> adt;
	std::vector<std::array<double, 4>> res;
	size_t size;

	CellPool(size_t num_cells)
		: pcell(num_cells)
		, pcell_new(num_cells)
		, p_q(num_cells)
		, qold(num_cells)
		, adt(num_cells)
		, res(num_cells)
		, size(num_cells)
	{}
};

struct NodePool {
	std::vector<std::array<double, 2>> p_x; // x
	size_t size;

	NodePool(size_t num_nodes)
		: p_x(num_nodes)
		, size(num_nodes)
	{}
};

inline void adt_calc(CellPool& cells) {
	for (size_t i = 0; i < cells.size; i++) {
		const double *x1 = cells.pcell_new[i][0].data();
		const double *x2 = cells.pcell_new[i][1].data();
		const double *x3 = cells.pcell_new[i][2].data();
		const double *x4 = cells.pcell_new[i][3].data();
		const double *q = cells.p_q[i].data();
		double *adt = &cells.adt[i];
		double dx, dy, ri, u, v, c;

		ri = 1.0f / q[0];
		u = ri * q[1];
		v = ri * q[2];
		c = sqrt(gam * gm1 * (ri * q[3] - 0.5f * (u * u + v * v)));

		dx = x2[0] - x1[0];
		dy = x2[1] - x1[1];
		*adt = fabs(u * dy - v * dx) + c * sqrt(dx * dx + dy * dy);

		dx = x3[0] - x2[0];
		dy = x3[1] - x2[1];
		*adt += fabs(u * dy - v * dx) + c * sqrt(dx * dx + dy * dy);

		dx = x4[0] - x3[0];
		dy = x4[1] - x3[1];
		*adt += fabs(u * dy - v * dx) + c * sqrt(dx * dx + dy * dy);

		dx = x1[0] - x4[0];
		dy = x1[1] - x4[1];
		*adt += fabs(u * dy - v * dx) + c * sqrt(dx * dx + dy * dy);

		*adt = (*adt) / cfl;
	}
}

inline void res_calc(const EdgePool& edges, const NodePool& nodes, CellPool& cells) {
	for (size_t i = 0; i < edges.size; i++) {
		const double *x1 = nodes.p_x[edges.pedge[i][0]].data();
		const double *x2 = nodes.p_x[edges.pedge[i][1]].data();
		const double *q1 = cells.p_q[edges.pecell[i][0]].data();
		const double *q2 = cells.p_q[edges.pecell[i][1]].data();
		const double *adt1 = &cells.adt[edges.pecell[i][0]];
		const double *adt2 = &cells.adt[edges.pecell[i][1]];
		double *res1 = cells.res[edges.pecell[i][0]].data();
		double *res2 = cells.res[edges.pecell[i][1]].data();

		double dx, dy, mu, ri, p1, vol1, p2, vol2, f;

		dx = x1[0] - x2[0];
		dy = x1[1] - x2[1];

		ri = 1.0f / q1[0];
		p1 = gm1 * (q1[3] - 0.5f * ri * (q1[1] * q1[1] + q1[2] * q1[2]));
		vol1 = ri * (q1[1] * dy - q1[2] * dx);

		ri = 1.0f / q2[0];
		p2 = gm1 * (q2[3] - 0.5f * ri * (q2[1] * q2[1] + q2[2] * q2[2]));
		vol2 = ri * (q2[1] * dy - q2[2] * dx);

		mu = 0.5f * ((*adt1) + (*adt2)) * eps;

		f = 0.5f * (vol1 * q1[0] + vol2 * q2[0]) + mu * (q1[0] - q2[0]);
		res1[0] += f;
		res2[0] -= f;
		f = 0.5f * (vol1 * q1[1] + p1 * dy + vol2 * q2[1] + p2 * dy) + mu * (q1[1] - q2[1]);
		res1[1] += f;
		res2[1] -= f;
		f = 0.5f * (vol1 * q1[2] - p1 * dx + vol2 * q2[2] - p2 * dx) + mu * (q1[2] - q2[2]);
		res1[2] += f;
		res2[2] -= f;
		f = 0.5f * (vol1 * (q1[3] + p1) + vol2 * (q2[3] + p2)) + mu * (q1[3] - q2[3]);
		res1[3] += f;
		res2[3] -= f;
	}
}

inline void bres_calc(const BackedgePool& bedges, const NodePool& nodes, CellPool& cells) {
	for (size_t i = 0; i < bedges.size; i++) {
		const double *x1 = nodes.p_x[bedges.pbedge[i][0]].data();
		const double *x2 = nodes.p_x[bedges.pbedge[i][1]].data();;
		const double *q1 = cells.p_q[bedges.pbecell[i]].data();
		const double *adt1 = &cells.adt[bedges.pbecell[i]];
		double *res1 = cells.res[bedges.pbecell[i]].data();
		const int *bound = &bedges.p_bound[i];

		double dx, dy, mu, ri, p1, vol1, p2, vol2, f;

		dx = x1[0] - x2[0];
		dy = x1[1] - x2[1];

		ri = 1.0f / q1[0];
		p1 = gm1 * (q1[3] - 0.5f * ri * (q1[1] * q1[1] + q1[2] * q1[2]));

		if (*bound == 1) {
			res1[1] += +p1 * dy;
			res1[2] += -p1 * dx;
		} else {
			vol1 = ri * (q1[1] * dy - q1[2] * dx);

			ri = 1.0f / qinf[0];
			p2 = gm1 * (qinf[3] - 0.5f * ri * (qinf[1] * qinf[1] + qinf[2] * qinf[2]));
			vol2 = ri * (qinf[1] * dy - qinf[2] * dx);

			mu = (*adt1) * eps;

			f = 0.5f * (vol1 * q1[0] + vol2 * qinf[0]) + mu * (q1[0] - qinf[0]);
			res1[0] += f;
			f = 0.5f * (vol1 * q1[1] + p1 * dy + vol2 * qinf[1] + p2 * dy) + mu * (q1[1] - qinf[1]);
			res1[1] += f;
			f = 0.5f * (vol1 * q1[2] - p1 * dx + vol2 * qinf[2] - p2 * dx) + mu * (q1[2] - qinf[2]);
			res1[2] += f;
			f = 0.5f * (vol1 * (q1[3] + p1) + vol2 * (qinf[3] + p2)) + mu * (q1[3] - qinf[3]);
			res1[3] += f;
		}
	}
}

inline void update(CellPool& cells, double* rms) {
	for (size_t i = 0; i < cells.size; i++) {
		const double *qold = cells.qold[i].data();
		double *q = cells.p_q[i].data();
		double *res = cells.res[i].data();
		const double *adt = &cells.adt[i];

		double del, adti;
		adti = 1.0f / (*adt);

		for (int n = 0; n < 4; n++) {
			del = adti * res[n];
			q[n] = qold[n] - del;
			res[n] = 0.0f;
			*rms += del * del;
		}
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

	size_t nnode, ncell, nedge, nbedge;
	if (fscanf(in_file, "%zu %zu %zu %zu \n", &nnode, &ncell, &nedge, &nbedge) != 4) {
		fprintf(stderr, "Could not parse size of lists of objects\n");
		return EXIT_FAILURE;
	}

	NodePool nodes(nnode);
	EdgePool edges(nedge);
	BackedgePool bedges(nbedge);
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

	for (size_t n = 0; n < nedge; n++) {
		if (fscanf(in_file, "%zu %zu %zu %zu \n",
				   &edges.pedge[n][0], &edges.pedge[n][1],
				   &edges.pecell[n][0], &edges.pecell[n][1]) != 4) {
			fprintf(stderr, "Error reading pedges and pecells\n");
			return EXIT_FAILURE;
		}
	}

	for (size_t n = 0; n < nbedge; n++) {
		if (fscanf(in_file, "%zu %zu %zu %d \n", &bedges.pbedge[n][0], &bedges.pbedge[n][1],
				   &bedges.pbecell[n], &bedges.p_bound[n]) != 4) {
			fprintf(stderr, "Error reading bedges, becells and bounds\n");
			return EXIT_FAILURE;
		}
	}
	fclose(in_file);

	for (size_t n = 0; n < ncell; n++) {
		for (size_t m = 0; m < 4; m++) {
			cells.p_q[n][m] = qinf[m];
			cells.res[n][m] = 0.0f;
		}
	}

	for (size_t i = 0; i < cells.size; i++) {
		for (size_t j = 0; j < 4; j++) {
			for (size_t k = 0; k < 4; k++) {
				cells.pcell_new[i][j][k] = nodes.p_x[cells.pcell[i][j]][k];
			}
		}
	}

	fprintf(stderr, "Parsing complete!\n");

	fprintf(stderr, "Now performing iterations...\n");

	struct timespec wall_start, wall_end, cpu_start, cpu_end;

	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpu_start);
	clock_gettime(CLOCK_REALTIME, &wall_start);

	double rms = 0;

	double adt_time = 0;
	double res_time = 0;
	double bres_time = 0;
	double upd_time = 0;

	for (int iter = 1; iter <= num_iter; iter++) {
		for (size_t i = 0; i < cells.size; i++) {
			cells.qold[i] = cells.p_q[i];
		}

		for (int k = 0; k < 2; k++) {

			struct timespec adt_start, adt_end;
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &adt_start);
			adt_calc(cells);
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &adt_end);
			adt_time += timespec_elapsed(&adt_end, &adt_start);

			struct timespec res_start, res_end;
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &res_start);
			res_calc(edges, nodes, cells);
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &res_end);
			res_time += timespec_elapsed(&res_end, &res_start);

			struct timespec bres_start, bres_end;
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &bres_start);
			bres_calc(bedges, nodes, cells);
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &bres_end);
			bres_time += timespec_elapsed(&bres_end, &bres_start);

			struct timespec upd_start, upd_end;
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &upd_start);
			update(cells, &rms);
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &upd_end);
			upd_time += timespec_elapsed(&upd_end, &upd_start);
		}

		rms = sqrt(rms / (double) cells.size);
	}

	printf("adt: %.8f, res: %.8f, bres: %.8f, upd: %.8f\n", adt_time, res_time, bres_time, upd_time);

	clock_gettime(CLOCK_REALTIME, &wall_end);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpu_end);

	fprintf(stderr, "Iterations complete!\n");

	double wall_elapsed = timespec_elapsed(&wall_end, &wall_start);
	printf("Wall clock time: %.8f\n", wall_elapsed);

	double cpu_elapsed = timespec_elapsed(&cpu_end, &cpu_start);
	printf("CPU clock time: %.8f\n", cpu_elapsed);

	fprintf(stderr, "Now writing to file %s...\n", argv[2]);

	FILE* out_file = fopen(argv[2], "w");
	for (size_t i = 0; i < cells.size; i++) {
		fprintf(out_file, "%f %f %f %f\n",
			cells.p_q[i][0], cells.p_q[i][1],
			cells.p_q[i][2], cells.p_q[i][3]);
	}
	fclose(out_file);

	return EXIT_SUCCESS;
}
