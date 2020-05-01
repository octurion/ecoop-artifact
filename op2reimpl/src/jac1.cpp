#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <array>
#include <vector>

static const size_t NN = 6;
static const int NUM_ITER = 2;

struct EdgePool {
	std::vector<std::array<size_t, 2>> ppedge;
	std::vector<std::array<double, 3>> p_A;
	size_t size;

	EdgePool(size_t num_edges)
		: ppedge(num_edges)
		, p_A(num_edges)
		, size(num_edges)
	{}
};

struct NodePool {
	std::vector<std::array<float, 2>> p_r;
	std::vector<std::array<float, 2>> p_u;
	std::vector<std::array<float, 3>> p_du;
	size_t size;

	NodePool(size_t num_nodes)
		: p_r(num_nodes)
		, p_u(num_nodes)
		, p_du(num_nodes)
		, size(num_nodes)
	{}
};

double timespec_elapsed(const struct timespec* end, const struct timespec* start)
{
	return end->tv_sec - start->tv_sec + 1e-9 * (end->tv_nsec - start->tv_nsec);
}

int main()
{
	size_t nnode = (NN - 1) * (NN - 1);
	size_t nedge = nnode + 4 * (NN - 1) * (NN - 2);

	NodePool nodes(nnode);
	EdgePool edges(nedge);

	fprintf(stderr, "Initialising dataset\n");

	size_t e = 0;

	for (size_t i = 1; i < NN; i++) {
		for (size_t j = 1; j < NN; j++) {
			size_t n = i - 1 + (j - 1) * (NN - 1);

			nodes.p_r[n][0] = 0;
			nodes.p_u[n][0] = 0;
			nodes.p_du[n][0] = 0;

			edges.ppedge[e][0] = n;
			edges.ppedge[e][1] = n;
			edges.p_A[e][0] = -1.0;
			e++;

			for (int pass = 0; pass < 4; pass++) {
				int i2 = i;
				int j2 = j;
				if (pass == 0) {
					i2 += -1;
				}
				if (pass == 1) {
					i2 += 1;
				}
				if (pass == 2) {
					j2 += -1;
				}
				if (pass == 3) {
					j2 += 1;
				}
				if ((i2 == 0) || (i2 == NN) || (j2 == 0) || (j2 == NN)) {
					nodes.p_r[n][0] += 0.25f;
				} else {
					edges.ppedge[e][0] = n;
					edges.ppedge[e][1] = i2 - 1 + (j2 - 1) * (NN - 1);
					edges.p_A[e][0] = 0.25;
					e++;
				}
			}
		}
	}

	fprintf(stderr, "Dataset complete!\n");

	fprintf(stderr, "Now performing iterations...\n");
	struct timespec wall_start, wall_end, cpu_start, cpu_end;

	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpu_start);
	clock_gettime(CLOCK_REALTIME, &wall_start);

	float alpha = 1.0f;

	float u_sum, u_max;
	float beta = 1.0;
	for (int iter = 0; iter < NUM_ITER; iter++) {
		for (size_t i = 0; i < edges.size; i++) {
			const double *A = &edges.p_A[i][0];
			const float *u = nodes.p_u[edges.ppedge[i][1]].data();
			float *du = nodes.p_du[edges.ppedge[i][0]].data();

			*du += (float)(beta * (*A) * (*u));
		}

		u_sum = 0;
		u_max = 0;
		for (size_t i = 0; i < nodes.size; i++) {
			const float *r = nodes.p_r[i].data();
			float *du = nodes.p_du[i].data();
			float *u = nodes.p_u[i].data();

			*u += *du + alpha * (*r);
			*du = 0.0f;
			u_sum += (*u) * (*u);
			u_max = std::max(u_max, *u);
		}
	}

	clock_gettime(CLOCK_REALTIME, &wall_end);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpu_end);
	fprintf(stderr, "Iterations complete!\n");

	for (size_t j = NN - 1; j > 0; j--) {
		for (size_t i = 1; i < NN; i++) {
			printf(" %7.4f", nodes.p_u[i - 1 + (j - 1) * (NN - 1)][0]);
		}
		printf("\n");
	}
	printf("\n");

	double wall_elapsed = timespec_elapsed(&wall_end, &wall_start);
	printf("Wall clock time: %.8f\n", wall_elapsed);

	double cpu_elapsed = timespec_elapsed(&cpu_end, &cpu_start);
	printf("CPU clock time: %.8f\n", cpu_elapsed);
	return EXIT_SUCCESS;
}
