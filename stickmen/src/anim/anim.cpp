#include "anim.h"

#include <xmmintrin.h>

static void animate_joints_scattered_owning_impl(JointScatteredOwning& joint)
{
	for (auto& e: joint.children) {
		e->orient_absolute = joint.orient_absolute * e->orient_relative;
		animate_joints_scattered_owning_impl(*e);
	}
}

void animate_joints_scattered_owning(JointScatteredOwning& root)
{
	root.orient_absolute = root.orient_relative;
	animate_joints_scattered_owning_impl(root);
}

static void animate_joints_scattered_owning_soa_impl(JointScatteredOwningSoa& joint)
{
	for (auto& e: joint.children) {
		e->orient_absolute = joint.orient_absolute * e->orient_relative;
		animate_joints_scattered_owning_soa_impl(*e);
	}
}

void animate_joints_scattered_owning_soa(JointScatteredOwningSoa& root)
{
	root.orient_absolute = root.orient_relative;
	animate_joints_scattered_owning_soa_impl(root);
}

static void animate_joints_scattered_ownerless_impl(JointScatteredOwnerless& joint)
{
	for (auto& e: joint.children) {
		e->orient_absolute = joint.orient_absolute * e->orient_relative;
		animate_joints_scattered_ownerless_impl(*e);
	}
}

void animate_joints_scattered_ownerless(JointScatteredOwnerless& root)
{
	root.orient_absolute = root.orient_relative;
	animate_joints_scattered_ownerless_impl(root);
}

static void animate_joints_scattered_ownerless_soa_impl(JointScatteredOwnerlessSoa& joint)
{
	for (auto& e: joint.children) {
		e->orient_absolute = joint.orient_absolute * e->orient_relative;
		animate_joints_scattered_ownerless_soa_impl(*e);
	}
}

void animate_joints_scattered_ownerless_soa(JointScatteredOwnerlessSoa& root)
{
	root.orient_absolute = root.orient_relative;
	animate_joints_scattered_ownerless_soa_impl(root);
}

void animate_joints_pooled_owning(std::vector<JointPooledOwning>& joints)
{
	joints[0].orient_relative = joints[0].orient_absolute;
	for (size_t i = 1; i < joints.size(); i++) {
		auto& e = joints[i];
		e.orient_absolute = e.parent->orient_relative * e.orient_relative;
	}
}

void animate_joints_pooled_owning_soa(std::vector<JointPooledOwningSoa>& joints)
{
	joints[0].orient_relative = joints[0].orient_absolute;
	for (size_t i = 1; i < joints.size(); i++) {
		auto& e = joints[i];
		e.orient_absolute = e.parent->orient_relative * e.orient_relative;
	}
}

void animate_joints_pooled_ownerless(std::vector<JointPooledOwnerless>& joints)
{
	joints[0].orient_relative = joints[0].orient_absolute;
	for (size_t i = 1; i < joints.size(); i++) {
		auto& e = joints[i];
		e.orient_absolute = e.parent->orient_relative * e.orient_relative;
	}
}

void animate_joints_pooled_ownerless_soa(std::vector<JointPooledOwnerlessSoa>& joints)
{
	joints[0].orient_relative = joints[0].orient_absolute;
	for (size_t i = 1; i < joints.size(); i++) {
		auto& e = joints[i];
		e.orient_absolute = e.parent->orient_relative * e.orient_relative;
	}
}

void animate_weights_scattered_parent(std::vector<WeightScatteredParent>& weights)
{
	for (auto& e: weights) {
		e.pos = e.joint->orient_absolute * e.pos;
	}
}

void animate_weights_pooled_parent(std::vector<WeightPooledParent>& weights)
{
	for (auto& e: weights) {
		e.pos = e.joint->orient_absolute * e.pos;
	}
}

void animate_weights_scattered_parent_soa(WeightScatteredParentSoa& weights)
{
	__m128 ONES = _mm_set1_ps(1);

	size_t len = weights.w.size();
	size_t cutoff = len & -4;
	for (size_t i = 0; i < cutoff; i += 4) {
		auto* joint0 = weights.joint[i + 0];
		auto* joint1 = weights.joint[i + 1];
		auto* joint2 = weights.joint[i + 2];
		auto* joint3 = weights.joint[i + 3];

		__m128 x = _mm_set_ps(
			joint3->orient_absolute.x,
			joint2->orient_absolute.x,
			joint1->orient_absolute.x,
			joint0->orient_absolute.x);

		__m128 y = _mm_set_ps(
			joint3->orient_absolute.y,
			joint2->orient_absolute.y,
			joint1->orient_absolute.y,
			joint0->orient_absolute.y);

		__m128 z = _mm_set_ps(
			joint3->orient_absolute.z,
			joint2->orient_absolute.z,
			joint1->orient_absolute.z,
			joint0->orient_absolute.z);

		__m128 w = _mm_set_ps(
			joint3->orient_absolute.w,
			joint2->orient_absolute.w,
			joint1->orient_absolute.w,
			joint0->orient_absolute.w);

		__m128 px = _mm_loadu_ps(&weights.posx[i]);
		__m128 py = _mm_loadu_ps(&weights.posy[i]);
		__m128 pz = _mm_loadu_ps(&weights.posz[i]);

		__m128 x2 = _mm_add_ps(x, x);
		__m128 y2 = _mm_add_ps(y, y);
		__m128 z2 = _mm_add_ps(z, z);
		__m128 xx2 = _mm_mul_ps(x, x2);
		__m128 yy2 = _mm_mul_ps(y, y2);
		__m128 zz2 = _mm_mul_ps(z, z2);

		__m128 xy2 = _mm_mul_ps(x, y2);
		__m128 wz2 = _mm_mul_ps(w, z2);
		__m128 xz2 = _mm_mul_ps(x, z2);
		__m128 wy2 = _mm_mul_ps(w, y2);
		__m128 yz2 = _mm_mul_ps(y, z2);
		__m128 wx2 = _mm_mul_ps(w, x2);

		__m128 a11 = _mm_sub_ps(_mm_sub_ps(ONES, yy2), zz2);
		__m128 a12 = _mm_add_ps(xy2, wz2);
		__m128 a13 = _mm_sub_ps(xz2, wy2);

		__m128 a21 = _mm_sub_ps(xy2, wz2);
		__m128 a22 = _mm_sub_ps(_mm_sub_ps(ONES, xx2), zz2);
		__m128 a23 = _mm_add_ps(yz2, wx2);

		__m128 a31 = _mm_add_ps(xz2, wy2);
		__m128 a32 = _mm_sub_ps(yz2, wx2);
		__m128 a33 = _mm_sub_ps(_mm_sub_ps(ONES, xx2), yy2);

		__m128 npx = _mm_mul_ps(a11, px);
		npx = _mm_add_ps(npx, _mm_mul_ps(a12, py));
		npx = _mm_add_ps(npx, _mm_mul_ps(a13, pz));

		__m128 npy = _mm_mul_ps(a21, py);
		npy = _mm_add_ps(npy, _mm_mul_ps(a22, py));
		npy = _mm_add_ps(npy, _mm_mul_ps(a23, pz));

		__m128 npz = _mm_mul_ps(a31, px);
		npz = _mm_add_ps(npz, _mm_mul_ps(a32, py));
		npz = _mm_add_ps(npz, _mm_mul_ps(a33, pz));

		_mm_storeu_ps(&weights.posx[i], npx);
		_mm_storeu_ps(&weights.posy[i], npy);
		_mm_storeu_ps(&weights.posz[i], npz);
	}

	for (size_t i = 0; i < len % 4; i++) {
		auto* joint = weights.joint[cutoff + i];

		vec3 pos;
		pos.x = weights.posx[cutoff + i];
		pos.y = weights.posy[cutoff + i];
		pos.z = weights.posz[cutoff + i];

		auto new_pos = joint->orient_absolute * pos;

		weights.posx[cutoff + i] = new_pos.x;
		weights.posy[cutoff + i] = new_pos.y;
		weights.posz[cutoff + i] = new_pos.z;
	}
}

void animate_weights_pooled_parent_soa(WeightPooledParentSoa& weights)
{
	__m128 ONES = _mm_set1_ps(1);

	size_t len = weights.w.size();
	size_t cutoff = len & -4;
	for (size_t i = 0; i < cutoff; i += 4) {
		auto* joint0 = weights.joint[i + 0];
		auto* joint1 = weights.joint[i + 1];
		auto* joint2 = weights.joint[i + 2];
		auto* joint3 = weights.joint[i + 3];

		__m128 x = _mm_set_ps(
			joint3->orient_absolute.x,
			joint2->orient_absolute.x,
			joint1->orient_absolute.x,
			joint0->orient_absolute.x);

		__m128 y = _mm_set_ps(
			joint3->orient_absolute.y,
			joint2->orient_absolute.y,
			joint1->orient_absolute.y,
			joint0->orient_absolute.y);

		__m128 z = _mm_set_ps(
			joint3->orient_absolute.z,
			joint2->orient_absolute.z,
			joint1->orient_absolute.z,
			joint0->orient_absolute.z);

		__m128 w = _mm_set_ps(
			joint3->orient_absolute.w,
			joint2->orient_absolute.w,
			joint1->orient_absolute.w,
			joint0->orient_absolute.w);

		__m128 px = _mm_loadu_ps(&weights.posx[i]);
		__m128 py = _mm_loadu_ps(&weights.posy[i]);
		__m128 pz = _mm_loadu_ps(&weights.posz[i]);

		__m128 x2 = _mm_add_ps(x, x);
		__m128 y2 = _mm_add_ps(y, y);
		__m128 z2 = _mm_add_ps(z, z);
		__m128 xx2 = _mm_mul_ps(x, x2);
		__m128 yy2 = _mm_mul_ps(y, y2);
		__m128 zz2 = _mm_mul_ps(z, z2);

		__m128 xy2 = _mm_mul_ps(x, y2);
		__m128 wz2 = _mm_mul_ps(w, z2);
		__m128 xz2 = _mm_mul_ps(x, z2);
		__m128 wy2 = _mm_mul_ps(w, y2);
		__m128 yz2 = _mm_mul_ps(y, z2);
		__m128 wx2 = _mm_mul_ps(w, x2);

		__m128 a11 = _mm_sub_ps(_mm_sub_ps(ONES, yy2), zz2);
		__m128 a12 = _mm_add_ps(xy2, wz2);
		__m128 a13 = _mm_sub_ps(xz2, wy2);

		__m128 a21 = _mm_sub_ps(xy2, wz2);
		__m128 a22 = _mm_sub_ps(_mm_sub_ps(ONES, xx2), zz2);
		__m128 a23 = _mm_add_ps(yz2, wx2);

		__m128 a31 = _mm_add_ps(xz2, wy2);
		__m128 a32 = _mm_sub_ps(yz2, wx2);
		__m128 a33 = _mm_sub_ps(_mm_sub_ps(ONES, xx2), yy2);

		__m128 npx = _mm_mul_ps(a11, px);
		npx = _mm_add_ps(npx, _mm_mul_ps(a12, py));
		npx = _mm_add_ps(npx, _mm_mul_ps(a13, pz));

		__m128 npy = _mm_mul_ps(a21, py);
		npy = _mm_add_ps(npy, _mm_mul_ps(a22, py));
		npy = _mm_add_ps(npy, _mm_mul_ps(a23, pz));

		__m128 npz = _mm_mul_ps(a31, px);
		npz = _mm_add_ps(npz, _mm_mul_ps(a32, py));
		npz = _mm_add_ps(npz, _mm_mul_ps(a33, pz));

		_mm_storeu_ps(&weights.posx[i], npx);
		_mm_storeu_ps(&weights.posy[i], npy);
		_mm_storeu_ps(&weights.posz[i], npz);
	}

	for (size_t i = 0; i < len % 4; i++) {
		auto* joint = weights.joint[cutoff + i];

		vec3 pos;
		pos.x = weights.posx[cutoff + i];
		pos.y = weights.posy[cutoff + i];
		pos.z = weights.posz[cutoff + i];

		auto new_pos = joint->orient_absolute * pos;

		weights.posx[cutoff + i] = new_pos.x;
		weights.posy[cutoff + i] = new_pos.y;
		weights.posz[cutoff + i] = new_pos.z;
	}
}

void animate_weights_scattered_parentless(JointScatteredOwning& root)
{
	for (auto& e: root.weights) {
		e.pos = root.orient_absolute * e.pos;
	}

	for (auto& e: root.children) {
		animate_weights_scattered_parentless(*e);
	}
}

void animate_weights_scattered_parentless_soa(JointScatteredOwningSoa& root)
{
	auto& q = root.orient_absolute;
	auto& weights = root.weights;

	float x = q.x;
	float y = q.y;
	float z = q.z;
	float w = q.w;

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

	__m128 a11 = _mm_set1_ps(1 - yy2 - zz2);
	__m128 a12 = _mm_set1_ps(xy2 + wz2);
	__m128 a13 = _mm_set1_ps(xz2 - wy2);

	__m128 a21 = _mm_set1_ps(xy2 - wz2);
	__m128 a22 = _mm_set1_ps(1 - xx2 - zz2);
	__m128 a23 = _mm_set1_ps(yz2 + wx2);

	__m128 a31 = _mm_set1_ps(xz2 + wy2);
	__m128 a32 = _mm_set1_ps(yz2 - wx2);
	__m128 a33 = _mm_set1_ps(1 - xx2 - yy2);

	size_t len = weights.w.size();
	size_t cutoff = len & -4;
	for (size_t i = 0; i < cutoff; i += 4) {
		__m128 px = _mm_loadu_ps(&weights.posx[i]);
		__m128 py = _mm_loadu_ps(&weights.posy[i]);
		__m128 pz = _mm_loadu_ps(&weights.posz[i]);

		__m128 npx = _mm_mul_ps(a11, px);
		npx = _mm_add_ps(npx, _mm_mul_ps(a12, py));
		npx = _mm_add_ps(npx, _mm_mul_ps(a13, pz));

		__m128 npy = _mm_mul_ps(a21, px);
		npy = _mm_add_ps(npy, _mm_mul_ps(a22, py));
		npy = _mm_add_ps(npy, _mm_mul_ps(a23, pz));

		__m128 npz = _mm_mul_ps(a31, px);
		npz = _mm_add_ps(npz, _mm_mul_ps(a32, py));
		npz = _mm_add_ps(npz, _mm_mul_ps(a33, pz));

		_mm_storeu_ps(&weights.posx[i], npx);
		_mm_storeu_ps(&weights.posy[i], npy);
		_mm_storeu_ps(&weights.posz[i], npz);
	}

	for (size_t i = 0; i < len % 4; i++) {
		vec3 pos;
		pos.x = weights.posx[cutoff + i];
		pos.y = weights.posy[cutoff + i];
		pos.z = weights.posz[cutoff + i];

		auto new_pos = q * pos;

		weights.posx[cutoff + i] = new_pos.x;
		weights.posy[cutoff + i] = new_pos.y;
		weights.posz[cutoff + i] = new_pos.z;
	}

	for (auto& e: root.children) {
		animate_weights_scattered_parentless_soa(*e);
	}
}

void animate_weights_pooled_parentless(std::vector<JointPooledOwning>& joints)
{
	for (auto& j: joints) {
		for (auto& e: j.weights) {
			e.pos = j.orient_absolute * e.pos;
		}
	}
}

void animate_weights_pooled_parentless_soa(std::vector<JointPooledOwningSoa>& joints)
{
	for (auto& j: joints) {
		auto& q = j.orient_absolute;
		auto& weights = j.weights;

		float x = q.x;
		float y = q.y;
		float z = q.z;
		float w = q.w;

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

		__m128 a11 = _mm_set1_ps(1 - yy2 - zz2);
		__m128 a12 = _mm_set1_ps(xy2 + wz2);
		__m128 a13 = _mm_set1_ps(xz2 - wy2);

		__m128 a21 = _mm_set1_ps(xy2 - wz2);
		__m128 a22 = _mm_set1_ps(1 - xx2 - zz2);
		__m128 a23 = _mm_set1_ps(yz2 + wx2);

		__m128 a31 = _mm_set1_ps(xz2 + wy2);
		__m128 a32 = _mm_set1_ps(yz2 - wx2);
		__m128 a33 = _mm_set1_ps(1 - xx2 - yy2);

		size_t len = weights.w.size();
		size_t cutoff = len & -4;
		for (size_t i = 0; i < cutoff; i += 4) {
			__m128 px = _mm_loadu_ps(&weights.posx[i]);
			__m128 py = _mm_loadu_ps(&weights.posy[i]);
			__m128 pz = _mm_loadu_ps(&weights.posz[i]);

			__m128 npx = _mm_mul_ps(a11, px);
			npx = _mm_add_ps(npx, _mm_mul_ps(a12, py));
			npx = _mm_add_ps(npx, _mm_mul_ps(a13, pz));

			__m128 npy = _mm_mul_ps(a21, px);
			npy = _mm_add_ps(npy, _mm_mul_ps(a22, py));
			npy = _mm_add_ps(npy, _mm_mul_ps(a23, pz));

			__m128 npz = _mm_mul_ps(a31, px);
			npz = _mm_add_ps(npz, _mm_mul_ps(a32, py));
			npz = _mm_add_ps(npz, _mm_mul_ps(a33, pz));

			_mm_storeu_ps(&weights.posx[i], npx);
			_mm_storeu_ps(&weights.posy[i], npy);
			_mm_storeu_ps(&weights.posz[i], npz);
		}

		for (size_t i = 0; i < len % 4; i++) {
			vec3 pos;
			pos.x = weights.posx[cutoff + i];
			pos.y = weights.posy[cutoff + i];
			pos.z = weights.posz[cutoff + i];

			auto new_pos = q * pos;

			weights.posx[cutoff + i] = new_pos.x;
			weights.posy[cutoff + i] = new_pos.y;
			weights.posz[cutoff + i] = new_pos.z;
		}
	}
}

void animate_joints_pooled_ownerless_scattered(std::vector<JointPooledOwnerlessScattered>& joints)
{
	joints[0].orient_relative = joints[0].orient_absolute;
	for (size_t i = 1; i < joints.size(); i++) {
		auto& e = joints[i];
		e.orient_absolute = e.parent->orient_relative * e.orient_relative;
	}
}

void animate_weights_pooled_parentless_scattered(std::vector<JointPooledOwningScattered>& joints)
{
	for (auto& j: joints) {
		for (auto& e: j.weights) {
			e->pos = j.orient_absolute * e->pos;
		}
	}
}

void animate_joints_pooled_owning_scattered(std::vector<JointPooledOwningScattered>& joints)
{
	joints[0].orient_relative = joints[0].orient_absolute;
	for (size_t i = 1; i < joints.size(); i++) {
		auto& e = joints[i];
		e.orient_absolute = e.parent->orient_relative * e.orient_relative;
	}
}

void animate_weights_pooled_parent_scattered(std::vector<std::unique_ptr<WeightPooledParentScattered>>& weights)
{
	for (auto& e: weights) {
		e->pos = e->joint->orient_absolute * e->pos;
	}
}
