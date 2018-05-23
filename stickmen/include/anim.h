#pragma once

#include <vector>

struct vec3
{
	float x, y, z;
};

struct quat
{
	float w, x, y, z;

	static inline quat from_xyz(float x, float y, float z)
	{
		quat r;

		r.x = x;
		r.y = y;
		r.z = z;
		r.w = 1 - x * x - y * y - z * z;

		return r;
	}

	inline friend quat operator*(const quat& a, const quat& b)
	{
		quat r;

		r.w = b.w * a.w - b.x * a.x - b.y * a.y - b.z * a.z;
		r.x = b.w * a.x + b.x * a.w - b.y * a.z + b.z * a.y;
		r.y = b.w * a.y + b.x * a.z + b.y * a.w - b.z * a.x;
		r.z = b.w * a.z + b.x * a.y + b.y * a.x - b.z * a.w;

		return r;
	}

	inline friend vec3 operator*(const quat& q, const vec3& v)
	{
		float x = q.x;
		float y = q.y;
		float z = q.z;
		float w = q.w;

		float px = v.x;
		float py = v.y;
		float pz = v.z;

		float a11 = 1 - 2*y*y - 2*z*z;
		float a12 = 2*x*y + 2*w*z;
		float a13 = 2*x*z - 2*w*y;

		float a21 = 2*x*y - 2*w*z;
		float a22 = 1 - 2*x*x - 2*z*z;
		float a23 = 2*y*z + 2*w*x;

		float a31 = 2*x*z + 2*w*y;
		float a32 = 2*y*z - 2*w*x;
		float a33 = 1 - 2*x*x - 2*y*y;

		vec3 r;

		r.x = a11*px + a12*py + a13*pz;
		r.y = a21*px + a22*py + a23*pz;
		r.z = a31*px + a32*py + a33*pz;

		return r;
	};
};

