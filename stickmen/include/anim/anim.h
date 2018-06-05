#pragma once

#include <memory>
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

		vec3 r;

		float x2 = x + x;
		float y2 = y + y;
		float z2 = z + z;
		float xx2 = x * x2;
		float yy2 = y * y2;
		float zz2 = z * z2;

		float xy2 = x * y2;
		float wz2 = w * z2;
		float xz2 = x * z2;
		float wy2 = w * y2;
		float yz2 = y * z2;
		float wx2 = w * x2;

		float a11 = 1 - yy2 - zz2;
		float a12 = xy2 + wz2;
		float a13 = xz2 + wy2;

		float a21 = xy2 - wz2;
		float a22 = 1 - xx2 - zz2;
		float a23 = yz2 + wx2;

		float a31 = xz2 - wy2;
		float a32 = yz2 - wx2;
		float a33 = 1 - xx2 - yy2;

		float npx = a11 * px + a12 * py + a13 * pz;
		float npy = a21 * px + a22 * py + a23 * pz;
		float npz = a31 * px + a32 * py + a33 * pz;

		r.x = npx;
		r.y = npy;
		r.z = npz;

		return r;
	};
};

struct JointScatteredOwning;
struct JointScatteredOwningSoa;
struct JointScatteredOwnerless;
struct JointScatteredOwnerlessSoa;

struct JointPooledOwning;
struct JointPooledOwningSoa;
struct JointPooledOwnerless;
struct JointPooledOwnerlessSoa;

struct JointPooledOwnerlessScattered;

struct WeightParentless;
struct WeightParentlessSoa;

struct WeightScatteredParent;
struct WeightScatteredParentSoa;

struct WeightPooledParent;
struct WeightPooledParentSoa;

struct WeightPooledParentScattered;

struct WeightParentless
{
	vec3 pos;
	float w;
};

struct WeightParentlessSoa
{
	std::vector<float> posx;
	std::vector<float> posy;
	std::vector<float> posz;

	std::vector<float> w;
};

struct WeightScatteredParent
{
	JointScatteredOwnerless* joint;

	vec3 pos;
	float w;
};

struct WeightScatteredParentSoa
{
	std::vector<JointScatteredOwnerlessSoa*> joint;

	std::vector<float> posx;
	std::vector<float> posy;
	std::vector<float> posz;

	std::vector<float> w;
};

struct WeightPooledParent
{
	JointPooledOwnerless* joint;

	vec3 pos;
	float w;
};

struct WeightPooledParentScattered
{
	JointPooledOwnerlessScattered* joint;

	vec3 pos;
	float w;
};

struct WeightPooledParentSoa
{
	std::vector<JointPooledOwnerlessSoa*> joint;

	std::vector<float> posx;
	std::vector<float> posy;
	std::vector<float> posz;

	std::vector<float> w;
};

struct JointScatteredOwning
{
	quat orient_relative;
	quat orient_absolute;

	std::vector<std::unique_ptr<JointScatteredOwning>> children;
	std::vector<WeightParentless> weights;
};

struct JointScatteredOwningSoa
{
	quat orient_relative;
	quat orient_absolute;

	std::vector<std::unique_ptr<JointScatteredOwningSoa>> children;
	WeightParentlessSoa weights;
};

struct JointScatteredOwnerless
{
	quat orient_relative;
	quat orient_absolute;

	std::vector<std::unique_ptr<JointScatteredOwnerless>> children;
};

struct JointScatteredOwnerlessSoa
{
	quat orient_relative;
	quat orient_absolute;

	std::vector<std::unique_ptr<JointScatteredOwnerlessSoa>> children;
};

struct JointPooledOwning
{
	JointPooledOwning* parent;

	quat orient_relative;
	quat orient_absolute;

	std::vector<WeightParentless> weights;
};

struct JointPooledOwningSoa
{
	JointPooledOwningSoa* parent;

	quat orient_relative;
	quat orient_absolute;

	WeightParentlessSoa weights;
};

struct JointPooledOwnerless
{
	JointPooledOwnerless* parent;

	quat orient_relative;
	quat orient_absolute;
};

struct JointPooledOwnerlessSoa
{
	JointPooledOwnerlessSoa* parent;

	quat orient_relative;
	quat orient_absolute;
};

struct JointPooledOwningScattered
{
	JointPooledOwningScattered* parent;

	quat orient_relative;
	quat orient_absolute;

	std::vector<std::unique_ptr<WeightParentless>> weights;
};

struct JointPooledOwnerlessScattered
{
	JointPooledOwnerlessScattered* parent;

	quat orient_relative;
	quat orient_absolute;
};

void animate_joints_scattered_owning(JointScatteredOwning& root);
void animate_joints_scattered_owning_soa(JointScatteredOwningSoa& root);
void animate_joints_scattered_ownerless(JointScatteredOwnerless& root);
void animate_joints_scattered_ownerless_soa(JointScatteredOwnerlessSoa& root);

void animate_joints_pooled_owning(std::vector<JointPooledOwning>& joints);
void animate_joints_pooled_owning_soa(std::vector<JointPooledOwningSoa>& joints);
void animate_joints_pooled_ownerless(std::vector<JointPooledOwnerless>& joints);
void animate_joints_pooled_ownerless_soa(std::vector<JointPooledOwnerlessSoa>& joints);

void animate_weights_scattered_parentless(JointScatteredOwning& root);
void animate_weights_scattered_parentless_soa(JointScatteredOwningSoa& root);

void animate_weights_scattered_parent(std::vector<WeightScatteredParent>& weights);
void animate_weights_scattered_parent_soa(WeightScatteredParentSoa& weights);

void animate_weights_pooled_parentless(std::vector<JointPooledOwning>& joints);
void animate_weights_pooled_parentless_soa(std::vector<JointPooledOwningSoa>& joints);

void animate_weights_pooled_parent(std::vector<WeightPooledParent>& weights);
void animate_weights_pooled_parent_soa(WeightPooledParentSoa& weights);

void animate_joints_pooled_ownerless_scattered(std::vector<JointPooledOwnerlessScattered>& joints);
void animate_weights_pooled_parentless_scattered(std::vector<JointPooledOwningScattered>& joints);

void animate_joints_pooled_owning_scattered(std::vector<JointPooledOwningScattered>& joints);
void animate_weights_pooled_parent_scattered(std::vector<std::unique_ptr<WeightPooledParentScattered>>& weights);
