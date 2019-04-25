#pragma once

#include <memory>
#include <vector>

struct vec3
{
	float x = 0, y = 0, z = 0;

	vec3 operator+(const vec3& rhs) const;
	vec3 operator-(const vec3& rhs) const;
};

struct quat
{
	float w = 1, x = 0, y = 0, z = 0;

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

struct JointMixed
{
	JointMixed* parent;
	JointMixed* next;
	vec3 pos;
	quat orient;

	std::vector<vec3> weight_initial_pos;
	std::vector<vec3> weight_pos;

	std::vector<float> weight_bias;

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

struct WeightPoolAos {
	std::vector<WeightAos> weights;
};

struct WeightPoolMixed {
	std::vector<vec3> initial_pos;
	std::vector<vec3> pos;
	std::vector<float> bias;
};

struct WeightPoolSoa {
	std::vector<float> initial_posx;
	std::vector<float> initial_posy;
	std::vector<float> initial_posz;

	std::vector<float> posx;
	std::vector<float> posy;
	std::vector<float> posz;

	std::vector<float> bias;
};

struct JointOnePool
{
	JointOnePool* parent;
	JointOnePool* next;

	vec3 pos;
	quat orient;

	size_t weights_start;
	size_t weights_end;

	void animate_my_weights(WeightPoolAos& weights);
	void animate_my_weights(WeightPoolMixed& weights);
	void animate_my_weights(WeightPoolSoa& weights);
};

struct WeightAosDynamic
{
	vec3 initial_pos;
	vec3 pos;
	float bias;

	WeightAosDynamic* next;
};

struct JointAosDynamic
{
	JointAosDynamic* parent;
	JointAosDynamic* next;

	vec3 pos;
	quat orient;

	WeightAosDynamic* weight_start;

	void animate_my_weights();
};

struct JointRaw
{
	size_t parent;
	size_t base_idx;

	vec3 base_pos;
	quat base_quat;
	int8_t flags;

	size_t idx_start;
	size_t num_weights;
	size_t idx_end;
};

struct WeightRaw
{
	size_t joint;
	vec3 pos;
	float bias;
};

struct Frame
{
	std::vector<float> values;
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
void animate_joints<JointMixed>(JointMixed* root, const JointRaw* base, const float* frame_components);

extern template
void animate_joints<JointSoa>(JointSoa* root, const JointRaw* base, const float* frame_components);

extern template
void animate_joints<JointOnePool>(JointOnePool* root, const JointRaw* base, const float* frame_components);

extern template
void animate_joints<JointAosDynamic>(JointAosDynamic* root, const JointRaw* base, const float* frame_components);

template<typename Joint>
void animate_weights(Joint* root);

extern template
void animate_weights<JointAos>(JointAos* root);

extern template
void animate_weights<JointMixed>(JointMixed* root);

extern template
void animate_weights<JointSoa>(JointSoa* root);

extern void animate_weights(JointAosDynamic* root);

void animate_weights(JointOnePool* root, WeightPoolAos& weights);
void animate_weights(JointOnePool* root, WeightPoolMixed& weights);
void animate_weights(JointOnePool* root, WeightPoolSoa& weights);
