#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <algorithm>
#include <array>
#include <unordered_set>
#include <vector>

static const size_t CHUNK_SIZE = 256;

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

struct Color {
	std::vector<std::pair<size_t, size_t>> chunk_ranges;
	std::unordered_set<size_t> indices;
};

struct Chunk {
	std::pair<size_t, size_t> range;
	std::vector<size_t> indices;
};

struct Edge {
	std::array<size_t, 2> pedge;  // Node indices
	std::array<size_t, 2> pecell; // Cell indices
};

struct Backedge {
	std::array<size_t, 2> pbedge; // Node indices
	size_t pbecell;               // Cell indices

	int p_bound; // Bound
};

struct Cell {
	std::array<size_t, 4> pcell; // Node indices
	// std::array<std::array<double, 2>, 4> pcell_new; // Node indices

	std::array<double, 4> p_q;
	std::array<double, 4> qold;
	double adt;
	std::array<double, 4> res;
};

struct Node {
	std::array<double, 2> p_x; // x
};

void copy_oldq(std::vector<Cell>& cells) {
	#pragma omp parallel for
	for (size_t i = 0; i < cells.size(); i++) {
		cells[i].qold = cells[i].p_q;
	}
}

void adt_calc(std::vector<Cell>& cells, const std::vector<Node>& nodes) {
	#pragma omp parallel for
	for (size_t i = 0; i < cells.size(); i++) {
		const double *x1 = nodes[cells[i].pcell[0]].p_x.data();
		const double *x2 = nodes[cells[i].pcell[1]].p_x.data();;
		const double *x3 = nodes[cells[i].pcell[2]].p_x.data();;
		const double *x4 = nodes[cells[i].pcell[3]].p_x.data();;
		const double *q = cells[i].p_q.data();
		double *adt = &cells[i].adt;
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

void res_calc(std::vector<Edge>& edges, const std::vector<Node>& nodes, std::vector<Cell>& cells, const std::vector<std::pair<size_t, size_t>>& ranges) {
	#pragma omp parallel
	for (const auto& r: ranges) {
		size_t begin = r.first;
		size_t end = r.second;

		#pragma omp for
		for (size_t i = begin; i < end; i++) {
			const double *x1 = nodes[edges[i].pedge[0]].p_x.data();
			const double *x2 = nodes[edges[i].pedge[1]].p_x.data();
			const double *q1 = cells[edges[i].pecell[0]].p_q.data();
			const double *q2 = cells[edges[i].pecell[1]].p_q.data();
			const double *adt1 = &cells[edges[i].pecell[0]].adt;
			const double *adt2 = &cells[edges[i].pecell[1]].adt;
			double *res1 = cells[edges[i].pecell[0]].res.data();
			double *res2 = cells[edges[i].pecell[1]].res.data();

			double q1_0 = q1[0];
			double q1_1 = q1[1];
			double q1_2 = q1[2];
			double q1_3 = q1[3];

			double q2_0 = q2[0];
			double q2_1 = q2[1];
			double q2_2 = q2[2];
			double q2_3 = q2[3];

			double ri1 = 1.0f / q1_0;
			double ri2 = 1.0f / q2_0;

			double dx = x1[0] - x2[0];
			double dy = x1[1] - x2[1];

			double p1 = gm1 * (q1_3 - 0.5f * ri1 * (q1_1 * q1_1 + q1_2 * q1_2));
			double vol1 = ri1 * (q1_1 * dy - q1_2 * dx);

			double p2 = gm1 * (q2_3 - 0.5f * ri2 * (q2_1 * q2_1 + q2_2 * q2_2));
			double vol2 = ri2 * (q2_1 * dy - q2_2 * dx);

			double mu = 0.5f * ((*adt1) + (*adt2)) * eps;

			double f;

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
}

void bres_calc(std::vector<Backedge>& bedges, const std::vector<Node>& nodes, std::vector<Cell>& cells) {
	for (size_t i = 0; i < bedges.size(); i++) {
		const double *x1 = nodes[bedges[i].pbedge[0]].p_x.data();
		const double *x2 = nodes[bedges[i].pbedge[1]].p_x.data();
		const double *q1 = cells[bedges[i].pbecell].p_q.data();
		const double *adt1 = &cells[bedges[i].pbecell].adt;
		const int *bound = &bedges[i].p_bound;
		double *res1 = cells[bedges[i].pbecell].res.data();

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

void update(std::vector<Cell>& cells, double* rms) {
	double rms_acc = 0;
	#pragma omp parallel for reduction(+:rms_acc)
	for (size_t i = 0; i < cells.size(); i++) {
		const double *qold = cells[i].qold.data();
		double *q = cells[i].p_q.data();
		double *res = cells[i].res.data();
		const double *adt = &cells[i].adt;

		double del, adti;
		adti = 1.0f / (*adt);

		for (int n = 0; n < 4; n++) {
			del = adti * res[n];
			q[n] = qold[n] - del;
			res[n] = 0.0f;
			rms_acc += del * del;
		}
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

	size_t nnode, ncell, nedge, nbedge;
	if (fscanf(in_file, "%zu %zu %zu %zu \n", &nnode, &ncell, &nedge, &nbedge) != 4) {
		fprintf(stderr, "Could not parse size of lists of objects\n");
		return EXIT_FAILURE;
	}

	std::vector<Node> nodes(nnode);
	std::vector<Edge> edges(nedge);
	std::vector<Backedge> bedges(nbedge);
	std::vector<Cell> cells(ncell);

	for (size_t n = 0; n < nnode; n++) {
		if (fscanf(in_file, "%lf %lf \n", &nodes[n].p_x[0], &nodes[n].p_x[1]) != 2) {
			fprintf(stderr, "Error reading nodes\n");
			return EXIT_FAILURE;
		}
	}

	for (size_t n = 0; n < ncell; n++) {
		if (fscanf(in_file, "%zu %zu %zu %zu \n",
				   &cells[n].pcell[0],
				   &cells[n].pcell[1],
				   &cells[n].pcell[2],
				   &cells[n].pcell[3]) != 4) {
			fprintf(stderr, "Error reading pcells\n");
			return EXIT_FAILURE;
		}
	}

	std::vector<Edge> edges_raw(nedge);

	for (size_t n = 0; n < nedge; n++) {
		if (fscanf(in_file, "%zu %zu %zu %zu \n",
				   &edges_raw[n].pedge[0], &edges_raw[n].pedge[1],
				   &edges_raw[n].pecell[0], &edges_raw[n].pecell[1]) != 4) {
			fprintf(stderr, "Error reading pedges and pecells\n");
			return EXIT_FAILURE;
		}
	}

	for (size_t n = 0; n < nbedge; n++) {
		if (fscanf(in_file, "%zu %zu %zu %d \n",
				   &bedges[n].pbedge[0], &bedges[n].pbedge[1],
				   &bedges[n].pbecell, &bedges[n].p_bound) != 4) {
			fprintf(stderr, "Error reading bedges, becells and bounds\n");
			return EXIT_FAILURE;
		}
	}
	fclose(in_file);

	std::vector<std::pair<size_t, size_t>> ranges;
	{
		std::vector<Color> colors;
		{
			size_t idx = 0;
			for (size_t chunk_id = 0; chunk_id < (nedge + CHUNK_SIZE - 1) / CHUNK_SIZE; chunk_id++) {
				Chunk chunk;
				chunk.range.first = idx;
				chunk.indices.reserve(2 * CHUNK_SIZE);
				for (size_t block_idx = 0; idx < nedge && block_idx < CHUNK_SIZE; idx++, block_idx++) {
					for (size_t k = 0; k < 2; k++) {
						chunk.indices.push_back(edges_raw[idx].pecell[k]);
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

		size_t last_idx = 0;
		size_t i = 0;
		for (const auto& c: colors) {
			for (const auto& r: c.chunk_ranges) {
				for (size_t k = r.first; k < r.second; k++, i++) {
					edges[i].pedge[0] = edges_raw[k].pedge[0];
					edges[i].pedge[1] = edges_raw[k].pedge[1];

					edges[i].pecell[0] = edges_raw[k].pecell[0];
					edges[i].pecell[1] = edges_raw[k].pecell[1];
				}
			}

			ranges.emplace_back(last_idx, i);
			last_idx = i;
		}
	}

	for (size_t n = 0; n < ncell; n++) {
		for (size_t m = 0; m < 4; m++) {
			cells[n].p_q[m] = qinf[m];
			cells[n].res[m] = 0.0f;
		}
	}

	fprintf(stderr, "Parsing complete!\n");

	fprintf(stderr, "Now performing iterations...\n");

	struct timespec wall_start, wall_end, cpu_start, cpu_end;

	clock_gettime(CLOCK_THREAD_CPUTIME_ID, &cpu_start);
	clock_gettime(CLOCK_REALTIME, &wall_start);

	double rms = 0;

	double adt_time = 0;
	double res_time = 0;
	double bres_time = 0;
	double upd_time = 0;

	for (int iter = 1; iter <= num_iter; iter++) {
		copy_oldq(cells);
		for (int k = 0; k < 2; k++) {
			adt_calc(cells, nodes);
			res_calc(edges, nodes, cells, ranges);
			bres_calc(bedges, nodes, cells);

			update(cells, &rms);
		}

		rms = sqrt(rms / (double) cells.size());

		if (iter % 100 == 0) {
			printf("iter = %4d, rms = %10.5e\n", iter, rms);
		}
	}

	printf("adt: %.8f, res: %.8f, bres: %.8f, upd: %.8f\n", adt_time, res_time, bres_time, upd_time);

	clock_gettime(CLOCK_REALTIME, &wall_end);
	clock_gettime(CLOCK_THREAD_CPUTIME_ID, &cpu_end);

	fprintf(stderr, "Iterations complete!\n");

	double wall_elapsed = timespec_elapsed(&wall_end, &wall_start);
	printf("Wall clock time: %.8f\n", wall_elapsed);

	double cpu_elapsed = timespec_elapsed(&cpu_end, &cpu_start);
	printf("CPU clock time: %.8f\n", cpu_elapsed);

	// fprintf(stderr, "Now writing to file %s...\n", argv[2]);

	// FILE* out_file = fopen(argv[2], "w");
	// for (size_t i = 0; i < cells.size(); i++) {
	// 	fprintf(out_file, "%f %f %f %f\n",
	// 		cells[i].p_q[0], cells[i].p_q[1],
	// 		cells[i].p_q[2], cells[i].p_q[3]);
	// }
	// fclose(out_file);

	return EXIT_SUCCESS;
}
