#include "anim.h"

#include <cmath>
#include <xmmintrin.h>

#define ANIM_XPOS (1 << 0)
#define ANIM_YPOS (1 << 1)
#define ANIM_ZPOS (1 << 2)

#define ANIM_XQUAT (1 << 3)
#define ANIM_YQUAT (1 << 4)
#define ANIM_ZQUAT (1 << 5)

template<typename Joint>
void animate_joints(Joint* root, const JointRaw* base, const float* frame_components)
{
	if (root == nullptr) {
		return;
	}

	size_t i = 0;
	for (Joint* it = root; it != nullptr; it = it->next, i++) {
		size_t idx = 0;

		auto flags = base[i].flags;
		auto base_pos = base[i].base_pos;
		auto base_orient = base[i].base_quat;
		auto base_idx = base[i].base_idx;

		if (flags & ANIM_XPOS) {
			base_pos.x = frame_components[base_idx + idx++];
		}
		if (flags & ANIM_YPOS) {
			base_pos.y = frame_components[base_idx + idx++];
		}
		if (flags & ANIM_ZPOS) {
			base_pos.z = frame_components[base_idx + idx++];
		}
		if (flags & ANIM_XQUAT) {
			base_orient.x = frame_components[base_idx + idx++];
		}
		if (flags & ANIM_YQUAT) {
			base_orient.y = frame_components[base_idx + idx++];
		}
		if (flags & ANIM_ZQUAT) {
			base_orient.z = frame_components[base_idx + idx++];
		}
		base_orient.gen_w();

		const auto* parent = it->parent;
		if (parent == nullptr) {
			it->pos = base_pos;
			it->orient = base_orient;
			continue;
		}

		auto rotated_pos = parent->orient.rotate(base_pos);
		it->pos = rotated_pos + parent->pos;
		it->orient = parent->orient * base_orient;
		it->orient.normalize();
	}
}

template<typename Joint>
void animate_weights(Joint* root)
{
	for (Joint* it = root; it != nullptr; it = it->next) {
		it->animate_my_weights();
	}
}

void animate_weights(JointOnePool* root, WeightPoolAos& weights)
{
	for (JointOnePool* it = root; it != nullptr; it = it->next) {
		it->animate_my_weights(weights);
	}
}
void animate_weights(JointOnePool* root, WeightPoolMixed& weights)
{
	for (JointOnePool* it = root; it != nullptr; it = it->next) {
		it->animate_my_weights(weights);
	}
}
void animate_weights(JointOnePool* root, WeightPoolSoa& weights)
{
	for (JointOnePool* it = root; it != nullptr; it = it->next) {
		it->animate_my_weights(weights);
	}
}

template
void animate_joints<JointAos>(JointAos* root, const JointRaw* base, const float* frame_components);

template
void animate_joints<JointMixed>(JointMixed* root, const JointRaw* base, const float* frame_components);

template
void animate_joints<JointSoa>(JointSoa* root, const JointRaw* base, const float* frame_components);

template
void animate_joints<JointOnePool>(JointOnePool* root, const JointRaw* base, const float* frame_components);

template
void animate_weights<JointAos>(JointAos* root);

template
void animate_weights<JointMixed>(JointMixed* root);

template
void animate_weights<JointSoa>(JointSoa* root);

inline void JointAos::animate_my_weights()
{
	auto num = weights.size();

	float x = orient.x;
	float y = orient.y;
	float z = orient.z;
	float w = orient.w;

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

	float posx = pos.x;
	float posy = pos.y;
	float posz = pos.z;

	for (size_t i = 0; i < num; i++) {
		float ix = weights[i].initial_pos.x;
		float iy = weights[i].initial_pos.y;
		float iz = weights[i].initial_pos.z;

		float npx = a11 * ix + a12 * iy + a13 * iz + posx;
		float npy = a21 * ix + a22 * iy + a23 * iz + posy;
		float npz = a31 * ix + a32 * iy + a33 * iz + posz;

		weights[i].pos.x = npx;
		weights[i].pos.y = npy;
		weights[i].pos.z = npz;
	}
}

inline void JointMixed::animate_my_weights()
{
	auto num = weight_initial_pos.size();

	float x = orient.x;
	float y = orient.y;
	float z = orient.z;
	float w = orient.w;

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

	float posx = pos.x;
	float posy = pos.y;
	float posz = pos.z;

	for (size_t i = 0; i < num; i++) {
		float ix = weight_initial_pos[i].x;
		float iy = weight_initial_pos[i].y;
		float iz = weight_initial_pos[i].z;

		float npx = a11 * ix + a12 * iy + a13 * iz + posx;
		float npy = a21 * ix + a22 * iy + a23 * iz + posy;
		float npz = a31 * ix + a32 * iy + a33 * iz + posz;

		weight_pos[i].x = npx;
		weight_pos[i].y = npy;
		weight_pos[i].z = npz;
	}
}

inline void JointSoa::animate_my_weights()
{
	auto num = weight_initial_posx.size();

	float x = orient.x;
	float y = orient.y;
	float z = orient.z;
	float w = orient.w;

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

	float posx = pos.x;
	float posy = pos.y;
	float posz = pos.z;

	for (size_t i = 0; i < num; i++) {
		float ix = weight_initial_posx[i];
		float iy = weight_initial_posy[i];
		float iz = weight_initial_posz[i];

		float npx = a11 * ix + a12 * iy + a13 * iz + posx;
		float npy = a21 * ix + a22 * iy + a23 * iz + posy;
		float npz = a31 * ix + a32 * iy + a33 * iz + posz;

		weight_posx[i] = npx;
		weight_posy[i] = npy;
		weight_posz[i] = npz;
	}
}

inline void JointOnePool::animate_my_weights(WeightPoolAos& wp)
{
	float x = orient.x;
	float y = orient.y;
	float z = orient.z;
	float w = orient.w;

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

	float posx = pos.x;
	float posy = pos.y;
	float posz = pos.z;

	for (size_t i = weights_start; i < weights_end; i++) {
		float ix = wp.weights[i].initial_pos.x;
		float iy = wp.weights[i].initial_pos.y;
		float iz = wp.weights[i].initial_pos.z;

		float npx = a11 * ix + a12 * iy + a13 * iz + posx;
		float npy = a21 * ix + a22 * iy + a23 * iz + posy;
		float npz = a31 * ix + a32 * iy + a33 * iz + posz;

		wp.weights[i].pos.x = npx;
		wp.weights[i].pos.y = npy;
		wp.weights[i].pos.z = npz;
	}
}

inline void JointOnePool::animate_my_weights(WeightPoolMixed& wp)
{
	float x = orient.x;
	float y = orient.y;
	float z = orient.z;
	float w = orient.w;

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

	float posx = pos.x;
	float posy = pos.y;
	float posz = pos.z;

	for (size_t i = weights_start; i < weights_end; i++) {
		float ix = wp.initial_pos[i].x;
		float iy = wp.initial_pos[i].y;
		float iz = wp.initial_pos[i].z;

		float npx = a11 * ix + a12 * iy + a13 * iz + posx;
		float npy = a21 * ix + a22 * iy + a23 * iz + posy;
		float npz = a31 * ix + a32 * iy + a33 * iz + posz;

		wp.pos[i].x = npx;
		wp.pos[i].y = npy;
		wp.pos[i].z = npz;
	}
}

inline void JointOnePool::animate_my_weights(WeightPoolSoa& wp)
{
	float x = orient.x;
	float y = orient.y;
	float z = orient.z;
	float w = orient.w;

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

	float posx = pos.x;
	float posy = pos.y;
	float posz = pos.z;

	for (size_t i = weights_start; i < weights_end; i++) {
		float ix = wp.initial_posx[i];
		float iy = wp.initial_posy[i];
		float iz = wp.initial_posz[i];

		float npx = a11 * ix + a12 * iy + a13 * iz + posx;
		float npy = a21 * ix + a22 * iy + a23 * iz + posy;
		float npz = a31 * ix + a32 * iy + a33 * iz + posz;

		wp.posx[i] = npx;
		wp.posy[i] = npy;
		wp.posz[i] = npz;
	}
}

inline vec3 vec3::operator+(const vec3& rhs) const
{
	vec3 ret;
	ret.x = x + rhs.x;
	ret.y = y + rhs.y;
	ret.z = z + rhs.z;

	return ret;
}

inline vec3 vec3::operator-(const vec3& rhs) const
{
	vec3 ret;
	ret.x = x - rhs.x;
	ret.y = y - rhs.y;
	ret.z = z - rhs.z;

	return ret;
}

void quat::gen_w()
{
	w = 1 - x * x - y * y - z * z;
	if (w < 0.0f) {
		w = 0.0f;
	} else {
		w = -sqrt(w);
	}
}

inline void quat::normalize()
{
	float mag = sqrt(x*x + y*y + z*z + w*w);
	if (mag < 1e-9f) {
		return;
	}

	float magRecip = 1.0f / mag;
	x *= magRecip;
	y *= magRecip;
	z *= magRecip;
	w *= magRecip;
}

inline quat operator*(const quat& a, const quat& b)
{
	quat r;

	r.w = b.w * a.w - b.x * a.x - b.y * a.y - b.z * a.z;
	r.x = b.w * a.x + b.x * a.w - b.y * a.z + b.z * a.y;
	r.y = b.w * a.y + b.x * a.z + b.y * a.w - b.z * a.x;
	r.z = b.w * a.z + b.x * a.y + b.y * a.x - b.z * a.w;

	return r;
}

inline vec3 quat::rotate(const vec3& v) const
{
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

	float npx = a11 * v.x + a12 * v.y + a13 * v.z;
	float npy = a21 * v.x + a22 * v.y + a23 * v.z;
	float npz = a31 * v.x + a32 * v.y + a33 * v.z;

	vec3 r;

	r.x = npx;
	r.y = npy;
	r.z = npz;

	return r;
}
