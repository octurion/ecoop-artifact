#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <emmintrin.h>

#include <array>
#include <vector>
#include <unordered_set>

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
	std::vector<size_t> pedge0;  // Node indices
	std::vector<size_t> pedge1;  // Node indices

	std::vector<size_t> pecell0; // Cell indices
	std::vector<size_t> pecell1; // Cell indices

	std::vector<double> x10;
	std::vector<double> x20;
	std::vector<double> x11;
	std::vector<double> x21;

	size_t size;

	EdgePool(size_t num_edges)
		: pedge0(num_edges)
		, pedge1(num_edges)
		, pecell0(num_edges)
		, pecell1(num_edges)
		, x10(num_edges)
		, x20(num_edges)
		, x11(num_edges)
		, x21(num_edges)
		, size(num_edges)
	{}
};

struct BackedgePool {
	std::vector<size_t> pbedge0; // Node indices
	std::vector<size_t> pbedge1; // Node indices
	std::vector<size_t> pbecell; // Cell indices

	std::vector<bool> p_bound; // Bound
	size_t size;

	BackedgePool(size_t num_backedges)
		: pbedge0(num_backedges)
		, pbedge1(num_backedges)
		, pbecell(num_backedges)
		, p_bound(num_backedges)
		, size(num_backedges)
	{}
};

struct CellPool {
	std::vector<size_t> pcell0; // Node indices
	std::vector<size_t> pcell1; // Node indices
	std::vector<size_t> pcell2; // Node indices
	std::vector<size_t> pcell3; // Node indices

	std::vector<double> p_x10;
	std::vector<double> p_x20;
	std::vector<double> p_x30;
	std::vector<double> p_x40;

	std::vector<double> p_x11;
	std::vector<double> p_x21;
	std::vector<double> p_x31;
	std::vector<double> p_x41;

	std::vector<double> p_q0;
	std::vector<double> p_q1;
	std::vector<double> p_q2;
	std::vector<double> p_q3;

	std::vector<double> qold0;
	std::vector<double> qold1;
	std::vector<double> qold2;
	std::vector<double> qold3;

	std::vector<double> adt;

	std::vector<double> res0;
	std::vector<double> res1;
	std::vector<double> res2;
	std::vector<double> res3;

	size_t size;

	CellPool(size_t num_cells)
		: pcell0(num_cells)
		, pcell1(num_cells)
		, pcell2(num_cells)
		, pcell3(num_cells)
		, p_x10(num_cells)
		, p_x20(num_cells)
		, p_x30(num_cells)
		, p_x40(num_cells)
		, p_x11(num_cells)
		, p_x21(num_cells)
		, p_x31(num_cells)
		, p_x41(num_cells)
		, p_q0(num_cells)
		, p_q1(num_cells)
		, p_q2(num_cells)
		, p_q3(num_cells)
		, qold0(num_cells)
		, qold1(num_cells)
		, qold2(num_cells)
		, qold3(num_cells)
		, adt(num_cells)
		, res0(num_cells)
		, res1(num_cells)
		, res2(num_cells)
		, res3(num_cells)
		, size(num_cells)
	{}
};

struct NodePool {
	std::vector<double> p_x0;
	std::vector<double> p_x1;
	size_t size;

	NodePool(size_t num_nodes)
		: p_x0(num_nodes)
		, p_x1(num_nodes)
		, size(num_nodes)
	{}
};

inline __m128d fabs(__m128d x) {
	__m128d sign_mask = _mm_set1_pd(-0.f);
	return _mm_andnot_pd(sign_mask, x);
}

inline void adt_calc(CellPool& cells) {
	for (size_t i = 0; i < (cells.size & -2); i += 2) {
		auto x10 = _mm_loadu_pd(&cells.p_x10[i]);
		auto x20 = _mm_loadu_pd(&cells.p_x20[i]);
		auto x30 = _mm_loadu_pd(&cells.p_x30[i]);
		auto x40 = _mm_loadu_pd(&cells.p_x40[i]);

		auto x11 = _mm_loadu_pd(&cells.p_x11[i]);
		auto x21 = _mm_loadu_pd(&cells.p_x21[i]);
		auto x31 = _mm_loadu_pd(&cells.p_x31[i]);
		auto x41 = _mm_loadu_pd(&cells.p_x41[i]);

		auto q0 = _mm_loadu_pd(&cells.p_q0[i]);
		auto q1 = _mm_loadu_pd(&cells.p_q1[i]);
		auto q2 = _mm_loadu_pd(&cells.p_q2[i]);
		auto q3 = _mm_loadu_pd(&cells.p_q3[i]);

		//double *adt = &cells.adt[i];
		//
		__m128d adt = _mm_setzero_pd();

		auto ri = _mm_set1_pd(1.0f) / q0;
		auto u = ri * q1;
		auto v = ri * q2;
		auto c = _mm_sqrt_pd(_mm_set1_pd(gam * gm1) * (ri * q3 - _mm_set1_pd(0.5f) * (u * u + v * v)));

		auto dx = x20 - x10;
		auto dy = x21 - x11;
		adt = fabs(u * dy - v * dx) + c * _mm_sqrt_pd(dx * dx + dy * dy);

		dx = x30 - x20;
		dy = x31 - x21;
		adt += fabs(u * dy - v * dx) + c * _mm_sqrt_pd(dx * dx + dy * dy);

		dx = x40 - x30;
		dy = x41 - x31;
		adt += fabs(u * dy - v * dx) + c * _mm_sqrt_pd(dx * dx + dy * dy);

		dx = x10 - x40;
		dy = x11 - x41;
		adt += fabs(u * dy - v * dx) + c * _mm_sqrt_pd(dx * dx + dy * dy);

		adt /= _mm_set1_pd(cfl);

		_mm_storeu_pd(&cells.adt[i], adt);
	}
	for (size_t i = (cells.size & -2); i < cells.size; i++) {
		const double x10 = cells.p_x10[i];
		const double x20 = cells.p_x20[i];
		const double x30 = cells.p_x30[i];
		const double x40 = cells.p_x40[i];

		const double x11 = cells.p_x11[i];
		const double x21 = cells.p_x21[i];
		const double x31 = cells.p_x31[i];
		const double x41 = cells.p_x41[i];

		const double q0 = cells.p_q0[i];
		const double q1 = cells.p_q1[i];
		const double q2 = cells.p_q2[i];
		const double q3 = cells.p_q3[i];

		double *adt = &cells.adt[i];
		double dx, dy, ri, u, v, c;

		ri = 1.0f / q0;
		u = ri * q1;
		v = ri * q2;
		c = sqrt(gam * gm1 * (ri * q3 - 0.5f * (u * u + v * v)));

		dx = x20 - x10;
		dy = x21 - x11;
		*adt = fabs(u * dy - v * dx) + c * sqrt(dx * dx + dy * dy);

		dx = x30 - x20;
		dy = x31 - x21;
		*adt += fabs(u * dy - v * dx) + c * sqrt(dx * dx + dy * dy);

		dx = x40 - x30;
		dy = x41 - x31;
		*adt += fabs(u * dy - v * dx) + c * sqrt(dx * dx + dy * dy);

		dx = x10 - x40;
		dy = x11 - x41;
		*adt += fabs(u * dy - v * dx) + c * sqrt(dx * dx + dy * dy);

		*adt = (*adt) / cfl;
	}
}

inline void res_calc(const EdgePool& edges, const NodePool& nodes, CellPool& cells) {
	for (size_t i = 0; i < (edges.size & -2); i += 2) {
		__m128d x10 = _mm_loadu_pd(&edges.x10[i]);
		__m128d x20 = _mm_loadu_pd(&edges.x20[i]);
		__m128d x11 = _mm_loadu_pd(&edges.x11[i]);
		__m128d x21 = _mm_loadu_pd(&edges.x21[i]);

		__m128d q10 = _mm_set_pd(
			cells.p_q0[edges.pecell0[i + 1]],
			cells.p_q0[edges.pecell0[i    ]]);
		__m128d q11 = _mm_set_pd(
			cells.p_q1[edges.pecell0[i + 1]],
			cells.p_q1[edges.pecell0[i    ]]);
		__m128d q12 = _mm_set_pd(
			cells.p_q2[edges.pecell0[i + 1]],
			cells.p_q2[edges.pecell0[i    ]]);
		__m128d q13 = _mm_set_pd(
			cells.p_q3[edges.pecell0[i + 1]],
			cells.p_q3[edges.pecell0[i    ]]);

		__m128d q20 = _mm_set_pd(
			cells.p_q0[edges.pecell1[i + 1]],
			cells.p_q0[edges.pecell1[i    ]]);
		__m128d q21 = _mm_set_pd(
			cells.p_q1[edges.pecell1[i + 1]],
			cells.p_q1[edges.pecell1[i    ]]);
		__m128d q22 = _mm_set_pd(
			cells.p_q2[edges.pecell1[i + 1]],
			cells.p_q2[edges.pecell1[i    ]]);
		__m128d q23 = _mm_set_pd(
			cells.p_q3[edges.pecell1[i + 1]],
			cells.p_q3[edges.pecell1[i    ]]);

		__m128d adt1 = _mm_set_pd(
			cells.adt[edges.pecell0[i + 1]],
			cells.adt[edges.pecell0[i   ]]);
		__m128d adt2 = _mm_set_pd(
			cells.adt[edges.pecell1[i + 1]],
			cells.adt[edges.pecell1[i   ]]);

		auto dx = x10 - x20;
		auto dy = x11 - x21;

		auto ri = _mm_set1_pd(1.0f) / q10;
		auto p1 = _mm_set1_pd(gm1) * (q13 - _mm_set1_pd(0.5f) * ri * (q11 * q11 + q12 * q12));
		auto vol1 = ri * (q11 * dy - q12 * dx);

		ri = _mm_set1_pd(1.0f) / q20;
		auto p2 = _mm_set1_pd(gm1) * (q23 - _mm_set1_pd(0.5f) * ri * (q21 * q21 + q22 * q22));
		auto vol2 = ri * (q21 * dy - q22 * dx);

		auto mu = _mm_set1_pd(0.5f) * (adt1 + adt2) * _mm_set1_pd(eps);

		auto f1 = _mm_set1_pd(0.5f) * (vol1 * q10 + vol2 * q20) + mu * (q10 - q20);
		auto f2 = _mm_set1_pd(0.5f) * (vol1 * q11 + p1 * dy + vol2 * q21 + p2 * dy) + mu * (q11 - q21);
		auto f3 = _mm_set1_pd(0.5f) * (vol1 * q12 - p1 * dx + vol2 * q22 - p2 * dx) + mu * (q12 - q22);
		auto f4 = _mm_set1_pd(0.5f) * (vol1 * (q13 + p1) + vol2 * (q23 + p2)) + mu * (q13 - q23);

		double *res10_1 = &cells.res0[edges.pecell0[i + 1]];
		double *res11_1 = &cells.res1[edges.pecell0[i + 1]];
		double *res12_1 = &cells.res2[edges.pecell0[i + 1]];
		double *res13_1 = &cells.res3[edges.pecell0[i + 1]];

		double *res10_0 = &cells.res0[edges.pecell0[i]];
		double *res11_0 = &cells.res1[edges.pecell0[i]];
		double *res12_0 = &cells.res2[edges.pecell0[i]];
		double *res13_0 = &cells.res3[edges.pecell0[i]];

		double *res20_1 = &cells.res0[edges.pecell1[i + 1]];
		double *res21_1 = &cells.res1[edges.pecell1[i + 1]];
		double *res22_1 = &cells.res2[edges.pecell1[i + 1]];
		double *res23_1 = &cells.res3[edges.pecell1[i + 1]];

		double *res20_0 = &cells.res0[edges.pecell1[i]];
		double *res21_0 = &cells.res1[edges.pecell1[i]];
		double *res22_0 = &cells.res2[edges.pecell1[i]];
		double *res23_0 = &cells.res3[edges.pecell1[i]];

		auto f1_0 = _mm_cvtsd_f64(f1);
		auto f1_1 = _mm_cvtsd_f64(_mm_shuffle_pd(f1, f1, 3));

		*res10_1 += f1_1;
		*res10_0 += f1_0;

		*res20_1 -= f1_1;
		*res20_0 -= f1_0;

		auto f2_0 = _mm_cvtsd_f64(f2);
		auto f2_1 = _mm_cvtsd_f64(_mm_shuffle_pd(f2, f2, 3));

		*res11_1 += f2_1;
		*res11_0 += f2_0;

		*res21_1 -= f2_1;
		*res21_0 -= f2_0;

		auto f3_0 = _mm_cvtsd_f64(f3);
		auto f3_1 = _mm_cvtsd_f64(_mm_shuffle_pd(f3, f3, 3));

		*res12_1 += f3_1;
		*res12_0 += f3_0;

		*res22_1 -= f3_1;
		*res22_0 -= f3_0;

		auto f4_0 = _mm_cvtsd_f64(f4);
		auto f4_1 = _mm_cvtsd_f64(_mm_shuffle_pd(f4, f4, 3));

		*res13_1 += f4_1;
		*res13_0 += f4_0;

		*res23_1 -= f4_1;
		*res23_0 -= f4_0;
	}
	for (size_t i = (edges.size & -2); i < edges.size; i++) {
		const double x10 = nodes.p_x0[edges.pedge0[i]];
		const double x20 = nodes.p_x0[edges.pedge1[i]];
		const double x11 = nodes.p_x1[edges.pedge0[i]];
		const double x21 = nodes.p_x1[edges.pedge1[i]];

		const double q10 = cells.p_q0[edges.pecell0[i]];
		const double q11 = cells.p_q1[edges.pecell0[i]];
		const double q12 = cells.p_q2[edges.pecell0[i]];
		const double q13 = cells.p_q3[edges.pecell0[i]];

		const double q20 = cells.p_q0[edges.pecell1[i]];
		const double q21 = cells.p_q1[edges.pecell1[i]];
		const double q22 = cells.p_q2[edges.pecell1[i]];
		const double q23 = cells.p_q3[edges.pecell1[i]];

		const double *adt1 = &cells.adt[edges.pecell0[i]];
		const double *adt2 = &cells.adt[edges.pecell1[i]];

		double *res10 = &cells.res0[edges.pecell0[i]];
		double *res11 = &cells.res1[edges.pecell0[i]];
		double *res12 = &cells.res2[edges.pecell0[i]];
		double *res13 = &cells.res3[edges.pecell0[i]];

		double *res20 = &cells.res0[edges.pecell1[i]];
		double *res21 = &cells.res1[edges.pecell1[i]];
		double *res22 = &cells.res2[edges.pecell1[i]];
		double *res23 = &cells.res3[edges.pecell1[i]];

		double dx, dy, mu, ri, p1, vol1, p2, vol2, f;

		dx = x10 - x20;
		dy = x11 - x21;

		ri = 1.0f / q10;
		p1 = gm1 * (q13 - 0.5f * ri * (q11 * q11 + q12 * q12));
		vol1 = ri * (q11 * dy - q12 * dx);

		ri = 1.0f / q20;
		p2 = gm1 * (q23 - 0.5f * ri * (q21 * q21 + q22 * q22));
		vol2 = ri * (q21 * dy - q22 * dx);

		mu = 0.5f * ((*adt1) + (*adt2)) * eps;

		f = 0.5f * (vol1 * q10 + vol2 * q20) + mu * (q10 - q20);
		*res10 += f;
		*res20 -= f;
		f = 0.5f * (vol1 * q11 + p1 * dy + vol2 * q21 + p2 * dy) + mu * (q11 - q21);
		*res11 += f;
		*res21 -= f;
		f = 0.5f * (vol1 * q12 - p1 * dx + vol2 * q22 - p2 * dx) + mu * (q12 - q22);
		*res12 += f;
		*res22 -= f;
		f = 0.5f * (vol1 * (q13 + p1) + vol2 * (q23 + p2)) + mu * (q13 - q23);
		*res13 += f;
		*res23 -= f;
	}
}

inline void bres_calc(const BackedgePool& bedges, const NodePool& nodes, CellPool& cells) {
	for (size_t i = 0; i < bedges.size; i++) {
		double x10 = nodes.p_x0[bedges.pbedge0[i]];
		double x11 = nodes.p_x1[bedges.pbedge0[i]];

		double x20 = nodes.p_x0[bedges.pbedge1[i]];
		double x21 = nodes.p_x1[bedges.pbedge1[i]];

		double q10 = cells.p_q0[bedges.pbecell[i]];
		double q11 = cells.p_q1[bedges.pbecell[i]];
		double q12 = cells.p_q2[bedges.pbecell[i]];
		double q13 = cells.p_q3[bedges.pbecell[i]];

		const double *adt1 = &cells.adt[bedges.pbecell[i]];

		double *res10 = &cells.res0[bedges.pbecell[i]];
		double *res11 = &cells.res1[bedges.pbecell[i]];
		double *res12 = &cells.res2[bedges.pbecell[i]];
		double *res13 = &cells.res3[bedges.pbecell[i]];

		double dx, dy, mu, ri, p1, vol1, p2, vol2, f;

		dx = x10 - x20;
		dy = x11 - x21;

		ri = 1.0f / q10;
		p1 = gm1 * (q13 - 0.5f * ri * (q11 * q11 + q12 * q12));

		if (bedges.p_bound[i]) {
			*res11 += +p1 * dy;
			*res12 += -p1 * dx;
		} else {
			vol1 = ri * (q11 * dy - q12 * dx);

			ri = 1.0f / qinf[0];
			p2 = gm1 * (qinf[3] - 0.5f * ri * (qinf[1] * qinf[1] + qinf[2] * qinf[2]));
			vol2 = ri * (qinf[1] * dy - qinf[2] * dx);

			mu = (*adt1) * eps;

			f = 0.5f * (vol1 * q10 + vol2 * qinf[0]) + mu * (q10 - qinf[0]);
			*res10 += f;
			f = 0.5f * (vol1 * q11 + p1 * dy + vol2 * qinf[1] + p2 * dy) + mu * (q11 - qinf[1]);
			*res11 += f;
			f = 0.5f * (vol1 * q12 - p1 * dx + vol2 * qinf[2] - p2 * dx) + mu * (q12 - qinf[2]);
			*res12 += f;
			f = 0.5f * (vol1 * (q13 + p1) + vol2 * (qinf[3] + p2)) + mu * (q13 - qinf[3]);
			*res13 += f;
		}
	}
}

inline void update(CellPool& cells, double* rms_ptr) {
	*rms_ptr = 0;
	__m128d rms = _mm_setzero_pd();
	for (size_t i = 0; i < (cells.size & -2); i += 2) {
		auto qold0 = _mm_loadu_pd(&cells.qold0[i]);
		double *q0 = &cells.p_q0[i];
		double *res0 = &cells.res0[i];

		auto adt = _mm_loadu_pd(&cells.adt[i]);
		auto adti = _mm_set1_pd(1.0) / adt;

		auto del = adti * _mm_loadu_pd(res0);
		_mm_storeu_pd(q0, qold0 - del);
		_mm_storeu_pd(res0, _mm_setzero_pd());
		rms += del * del;
	}
	for (size_t i = (cells.size & -2); i < cells.size; i++) {
		double qold0 = cells.qold0[i];
		double *q0 = &cells.p_q0[i];
		double *res0 = &cells.res0[i];

		double adt = cells.adt[i];
		double adti = 1.0 / adt;

		double del = adti * *res0;
		*q0 = qold0 - del;
		*res0 = 0.0f;
		*rms_ptr += del * del;
	}

	for (size_t i = 0; i < (cells.size & -2); i += 2) {
		auto qold1 = _mm_loadu_pd(&cells.qold1[i]);
		double *q1 = &cells.p_q1[i];
		double *res1 = &cells.res1[i];

		auto adt = _mm_loadu_pd(&cells.adt[i]);
		auto adti = _mm_set1_pd(1.0) / adt;

		auto del = adti * _mm_loadu_pd(res1);
		_mm_storeu_pd(q1, qold1 - del);
		_mm_storeu_pd(res1, _mm_setzero_pd());
		rms += del * del;
	}
	for (size_t i = (cells.size & -2); i < cells.size; i++) {
		double qold1 = cells.qold1[i];
		double *q1 = &cells.p_q1[i];
		double *res1 = &cells.res1[i];

		double adt = cells.adt[i];
		double adti = 1.0 / adt;

		double del = adti * *res1;
		*q1 = qold1 - del;
		*res1 = 0.0f;
		*rms_ptr += del * del;
	}

	for (size_t i = 0; i < (cells.size & -2); i += 2) {
		auto qold2 = _mm_loadu_pd(&cells.qold2[i]);
		double *q2 = &cells.p_q2[i];
		double *res2 = &cells.res2[i];

		auto adt = _mm_loadu_pd(&cells.adt[i]);
		auto adti = _mm_set1_pd(1.0) / adt;

		auto del = adti * _mm_loadu_pd(res2);
		_mm_storeu_pd(q2, qold2 - del);
		_mm_storeu_pd(res2, _mm_setzero_pd());
		rms += del * del;
	}
	for (size_t i = (cells.size & -2); i < cells.size; i++) {
		double qold2 = cells.qold2[i];
		double *q2 = &cells.p_q2[i];
		double *res2 = &cells.res2[i];

		double adt = cells.adt[i];
		double adti = 1.0 / adt;

		double del = adti * *res2;
		*q2 = qold2 - del;
		*res2 = 0.0f;
		*rms_ptr += del * del;
	}

	for (size_t i = 0; i < (cells.size & -2); i += 2) {
		auto qold3 = _mm_loadu_pd(&cells.qold3[i]);
		double *q3 = &cells.p_q3[i];
		double *res3 = &cells.res3[i];

		auto adt = _mm_loadu_pd(&cells.adt[i]);
		auto adti = _mm_set1_pd(1.0) / adt;

		auto del = adti * _mm_loadu_pd(res3);
		_mm_storeu_pd(q3, qold3 - del);
		_mm_storeu_pd(res3, _mm_setzero_pd());
		rms += del * del;
	}
	for (size_t i = (cells.size & -2); i < cells.size; i++) {
		double qold3 = cells.qold3[i];
		double *q3 = &cells.p_q3[i];
		double *res3 = &cells.res3[i];

		double adt = cells.adt[i];
		double adti = 1.0f / (adt);

		double del = adti * *res3;
		*q3 = qold3 - del;
		*res3 = 0.0f;
		*rms_ptr += del * del;
	}
	*rms_ptr += _mm_cvtsd_f64(rms);
	*rms_ptr += _mm_cvtsd_f64(_mm_shuffle_pd(rms, rms, 3));
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
		if (fscanf(in_file, "%lf %lf \n", &nodes.p_x0[n], &nodes.p_x1[n]) != 2) {
			fprintf(stderr, "Error reading nodes\n");
			return EXIT_FAILURE;
		}
	}

	for (size_t n = 0; n < ncell; n++) {
		if (fscanf(in_file, "%zu %zu %zu %zu \n",
				   &cells.pcell0[n],
				   &cells.pcell1[n],
				   &cells.pcell2[n],
				   &cells.pcell3[n]) != 4) {
			fprintf(stderr, "Error reading pcells\n");
			return EXIT_FAILURE;
		}
	}

	for (size_t n = 0; n < nedge; n++) {
		if (fscanf(in_file, "%zu %zu %zu %zu \n",
				   &edges.pedge0[n], &edges.pedge1[n],
				   &edges.pecell0[n], &edges.pecell1[n]) != 4) {
			fprintf(stderr, "Error reading pedges and pecells\n");
			return EXIT_FAILURE;
		}
	}

	for (size_t n = 0; n < nbedge; n++) {
		int bound;
		if (fscanf(in_file, "%zu %zu %zu %d \n", &bedges.pbedge0[n], &bedges.pbedge1[n],
				   &bedges.pbecell[n], &bound) != 4) {
			fprintf(stderr, "Error reading bedges, becells and bounds\n");
			return EXIT_FAILURE;
		}
		bedges.p_bound[n] = bound;
	}
	fclose(in_file);

	for (size_t i = 0; i < cells.size; i++) {
		cells.p_q0[i] = qinf[0];
		cells.p_q1[i] = qinf[1];
		cells.p_q2[i] = qinf[2];
		cells.p_q3[i] = qinf[3];

		cells.res0[i] = 0.0f;
		cells.res1[i] = 0.0f;
		cells.res2[i] = 0.0f;
		cells.res3[i] = 0.0f;

		cells.p_x10[i] = nodes.p_x0[cells.pcell0[i]];
		cells.p_x20[i] = nodes.p_x0[cells.pcell1[i]];
		cells.p_x30[i] = nodes.p_x0[cells.pcell2[i]];
		cells.p_x40[i] = nodes.p_x0[cells.pcell3[i]];

		cells.p_x11[i] = nodes.p_x1[cells.pcell0[i]];
		cells.p_x21[i] = nodes.p_x1[cells.pcell1[i]];
		cells.p_x31[i] = nodes.p_x1[cells.pcell2[i]];
		cells.p_x41[i] = nodes.p_x1[cells.pcell3[i]];
	}
	for (size_t i = 0; i < edges.size; i++) {
		edges.x10[i] = nodes.p_x0[edges.pedge0[i]];
		edges.x20[i] = nodes.p_x0[edges.pedge1[i]];
		edges.x11[i] = nodes.p_x1[edges.pedge0[i]];
		edges.x21[i] = nodes.p_x1[edges.pedge1[i]];
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
			cells.qold0[i] = cells.p_q0[i];
			cells.qold1[i] = cells.p_q1[i];
			cells.qold2[i] = cells.p_q2[i];
			cells.qold3[i] = cells.p_q3[i];
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
			cells.p_q0[i], cells.p_q1[i],
			cells.p_q2[i], cells.p_q3[i]);
	}
	fclose(out_file);

	return EXIT_SUCCESS;
}
