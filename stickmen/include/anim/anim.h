#pragma once

#include <memory>
#include <vector>

struct vec3
{
	float x, y, z;

	vec3 operator+(const vec3& rhs) const;
	vec3 operator-(const vec3& rhs) const;
};

struct quat
{
	float w, x, y, z;

	void gen_w();
	void normalize();
	vec3 rotate(const vec3& v) const;

	friend quat operator*(const quat& a, const quat& b);
};

struct WeightAos
{
	vec3 initial_pos;
	vec3 pos;
	float bias;
};

struct JointAos
{
	JointAos* parent;
	JointAos* next;

	vec3 pos;
	quat orient;

	std::vector<WeightAos> weights;

	void animate_my_weights();
};

struct JointSoa
{
	JointSoa* parent;
	JointSoa* next;
	vec3 pos;
	quat orient;

	std::vector<float> weight_initial_posx;
	std::vector<float> weight_initial_posy;
	std::vector<float> weight_initial_posz;

	std::vector<float> weight_posx;
	std::vector<float> weight_posy;
	std::vector<float> weight_posz;

	std::vector<float> weight_bias;

	void animate_my_weights();
};

struct JointRaw
{
	size_t parent;
	size_t base_idx;

	vec3 base_pos;
	quat base_quat;
	int8_t flags;
};

struct Frame
{
	std::vector<float> values;
};

struct WeightRaw
{
	size_t joint;
	vec3 pos;
	float bias;
};

struct ModelInput
{
	std::vector<JointRaw> joints_raw;
	std::vector<WeightRaw> weights_raw;

	std::vector<Frame> frames;
};

template<typename Joint>
void animate_joints(Joint* root, const JointRaw* base, const float* frame_components);

extern template
void animate_joints<JointAos>(JointAos* root, const JointRaw* base, const float* frame_components);

extern template
void animate_joints<JointSoa>(JointSoa* root, const JointRaw* base, const float* frame_components);

template<typename Joint>
void animate_weights(Joint* root);

extern template
void animate_weights<JointAos>(JointAos* root);

extern template
void animate_weights<JointSoa>(JointSoa* root);
