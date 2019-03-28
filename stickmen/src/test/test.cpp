#include "anim/anim.h"

#include "test/dataset.h"

#include <xmmintrin.h>

#include <benchmark/benchmark.h>

#include <deque>
#include <memory>
#include <numeric>
#include <utility>

#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#define NUM_MODELS 2
#define TO_REMOVE_RATIO 0.1f

struct BenchmarkOptions
{
	size_t num_copies;
	size_t stride;
	bool flush_range;
};

static BenchmarkOptions from_args(benchmark::State& state)
{
	BenchmarkOptions opt;

	opt.num_copies = state.range(0);
	opt.stride = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);
	opt.flush_range = state.range(1);

	return opt;
}

static bool initialized = false;
static ModelInput INPUTS[NUM_MODELS];

static void CustomArgs(benchmark::internal::Benchmark* b)
{
	for (int i = 1; i <= 10; i++) {
		b->Args({i * 500, 1});
	}
}

static void flush_range(const void* ptr, size_t stride, size_t len)
{
	size_t flush_count = (len + stride) & (stride - 1);

	for (size_t i = 0; i < flush_count; i += stride) {
		_mm_clflush((const char*) ptr + i * stride);
	}
}

ModelInput parse_files(const char* mesh_filename, const char* anim_filename)
{
	ModelInput input;

	size_t len = 0;
	char* line = NULL;

	FILE* mesh_in = fopen(mesh_filename, "r");
	if (mesh_in == NULL) {
		perror(mesh_filename);
		abort();
	}

	while (getline(&line, &len, mesh_in) != -1) {
		WeightRaw w;
		if (sscanf(line, " weight %*d %zu %f ( %f %f %f )",
				   &w.joint, &w.bias, &w.pos.x, &w.pos.y, &w.pos.z) != 5) {
			continue;
		}

		input.weights_raw.push_back(w);
	}
	fclose(mesh_in);

	FILE* anim_in = fopen(anim_filename, "r");
	if (anim_in == NULL) {
		perror(anim_filename);
		abort();
	}

	size_t num_anim_components;

	while (getline(&line, &len, anim_in) != -1) {
		if (sscanf(line, " numAnimatedComponents %zu", &num_anim_components) == 1) {
			break;
		}
	}

	while (getline(&line, &len, anim_in) != -1) {
		const char* pattern = "hierarchy {";
		if (strncmp(line, pattern, strlen(pattern)) == 0) {
			break;
		}
	}

	while (getline(&line, &len, anim_in) != -1) {
		const char* pattern = "}";
		if (strncmp(line, pattern, strlen(pattern)) == 0) {
			break;
		}

		JointRaw j;
		if (sscanf(line, " %*s %zu %" SCNi8 " %zu", &j.parent, &j.flags, &j.base_idx) != 3) {
			fprintf(stderr, "Couldn't read joint hierarchy\n");
			abort();
		}

		j.base_pos.x = 0;
		j.base_pos.y = 0;
		j.base_pos.z = 0;

		j.base_quat.x = 0;
		j.base_quat.y = 0;
		j.base_quat.z = 0;
		j.base_quat.w = 0;

		j.idx_start = 0;
		j.idx_end = 0;
		j.num_weights = 0;

		input.joints_raw.push_back(j);
	}

	while (getline(&line, &len, anim_in) != -1) {
		const char* pattern = "baseframe {";
		if (strncmp(line, pattern, strlen(pattern)) == 0) {
			break;
		}
	}

	size_t idx = 0;
	while (getline(&line, &len, anim_in) != -1) {
		const char* pattern = "}";
		if (strncmp(line, pattern, strlen(pattern)) == 0) {
			break;
		}

		auto& j = input.joints_raw[idx];
		if (sscanf(line, " ( %f %f %f ) ( %f %f %f )",
				   &j.base_pos.x, &j.base_pos.y, &j.base_pos.z,
				   &j.base_quat.x, &j.base_quat.y, &j.base_quat.z) != 6) {
			fprintf(stderr, "Couldn't read base frame\n");
			abort();
		}

		j.base_quat.gen_w();
		idx++;
	}

	if (idx != input.joints_raw.size()) {
		fprintf(stderr, "Couldn't read base frame\n");
		abort();
	}

	while (getline(&line, &len, anim_in) != -1) {
		int frame_id;
		if (sscanf(line, " frame %d {", &frame_id) != 1) {
			continue;
		}

		Frame frame;
		for (size_t i = 0; i < num_anim_components; i++) {
			float value;
			int ret = fscanf(anim_in, "%f", &value);
			(void) ret;

			frame.values.push_back(value);
		}

		input.frames.emplace_back(std::move(frame));
	}

	fclose(anim_in);

	std::stable_sort(
		input.weights_raw.begin(),
		input.weights_raw.end(),
		[](const WeightRaw& lhs, const WeightRaw& rhs) {
			return lhs.joint < rhs.joint;
		});

	for (const auto& e: input.weights_raw) {
		input.joints_raw[e.joint].num_weights++;
	}

	input.joints_raw[0].idx_end = input.joints_raw[0].num_weights;
	for (size_t i = 1; i < input.joints_raw.size(); i++) {
		input.joints_raw[i].idx_start = input.joints_raw[i - 1].idx_end;
		input.joints_raw[i].idx_end = input.joints_raw[i].idx_start +
			input.joints_raw[i].num_weights;
	}

	free(line);

	return input;
}

struct PooledJointsAosWeights
{
	int frame;
	const ModelInput* input;
	std::vector<JointAos> joints;
};

struct Bench_PooledJointsAosWeights
{
	std::deque<PooledJointsAosWeights> models;
	size_t counter = 0;

	PooledJointsAosWeights construct_model(const ModelInput& input, size_t id)
	{
		PooledJointsAosWeights model;
		model.input = &input;
		model.frame = id % input.frames.size();

		for (size_t i = 0; i < input.joints_raw.size(); i++) {
			JointAos j;
			j.orient = input.joints_raw[i].base_quat;
			j.pos = input.joints_raw[i].base_pos;

			model.joints.emplace_back(std::move(j));
		}

		model.joints[0].parent = nullptr;
		model.joints[0].next = nullptr;
		for (size_t i = 1; i < model.joints.size(); i++) {
			model.joints[i].parent = &model.joints[input.joints_raw[i].parent];
			model.joints[i].next = nullptr;
			model.joints[i - 1].next = &model.joints[i];
		}

		for (size_t i = 0; i < input.weights_raw.size(); i++) {
			WeightAos w;
			w.bias = input.weights_raw[i].bias;
			w.initial_pos = input.weights_raw[i].pos;
			w.pos.x = 0;
			w.pos.y = 0;
			w.pos.z = 0;

			model.joints[input.weights_raw[i].joint].weights.push_back(w);
		}

		return model;
	}

	void initialize_from_opts(BenchmarkOptions& opts)
	{
		for (size_t i = 0; i < opts.num_copies; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				models.emplace_back(construct_model(INPUTS[j], i));
			}
		}

		counter = opts.num_copies - 1;
	}

	void flush_from_cache(const BenchmarkOptions& opts)
	{
		for (const auto& m: models) {
			for (const auto& j: m.joints) {
				size_t range = j.weights.size() * sizeof(decltype(j.weights)::value_type);
				flush_range(j.weights.data(), opts.stride, range);
			}

			size_t range = m.joints.size() * sizeof(decltype(m.joints)::value_type);
			flush_range(m.joints.data(), opts.stride, range);
		}
	}

	void run_animate_joints()
	{
		for (auto& m: models) {
			animate_joints(&m.joints[0],
						   m.input->joints_raw.data(),
						   m.input->frames[m.frame].values.data());
		}
	}

	void run_animate_weights()
	{
		for (auto& m: models) {
			animate_weights(&m.joints[0]);

			m.frame++;
			m.frame %= m.input->frames.size();
		}
	}

	void modify_structure(const BenchmarkOptions& opts)
	{
		size_t num_to_replace = opts.num_copies * TO_REMOVE_RATIO;
		for (size_t i = 0; i < num_to_replace; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				models.pop_front();
			}
		}

		for (size_t i = 0; i < num_to_replace; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				models.emplace_back(construct_model(INPUTS[j], counter));
			}
			counter++;
		}
	}
};

struct ScatteredJointsAosWeights
{
	int frame;
	const ModelInput* input;
	std::vector<std::unique_ptr<JointAos>> joints;
};

struct Bench_ScatteredJointsAosWeights
{
	std::deque<ScatteredJointsAosWeights> models;
	size_t counter = 0;

	ScatteredJointsAosWeights construct_model(const ModelInput& input, size_t id)
	{
		ScatteredJointsAosWeights model;
		model.input = &input;
		model.frame = id % input.frames.size();

		for (size_t i = 0; i < input.joints_raw.size(); i++) {
			std::unique_ptr<JointAos> j(new JointAos);
			j->orient = input.joints_raw[i].base_quat;
			j->pos = input.joints_raw[i].base_pos;

			model.joints.emplace_back(std::move(j));
		}

		model.joints[0]->parent = nullptr;
		model.joints[0]->next = nullptr;
		for (size_t i = 1; i < model.joints.size(); i++) {
			model.joints[i]->parent = model.joints[input.joints_raw[i].parent].get();
			model.joints[i]->next = nullptr;
			model.joints[i - 1]->next = model.joints[i].get();
		}

		for (size_t i = 0; i < input.weights_raw.size(); i++) {
			WeightAos w;
			w.bias = input.weights_raw[i].bias;
			w.initial_pos = input.weights_raw[i].pos;
			w.pos.x = 0;
			w.pos.y = 0;
			w.pos.z = 0;

			model.joints[input.weights_raw[i].joint]->weights.push_back(w);
		}

		return model;
	}

	void initialize_from_opts(BenchmarkOptions& opts)
	{
		for (size_t i = 0; i < opts.num_copies; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				models.emplace_back(construct_model(INPUTS[j], i));
			}
		}

		counter = opts.num_copies - 1;
	}

	void flush_from_cache(const BenchmarkOptions& opts)
	{
		for (const auto& m: models) {
			for (const auto& j: m.joints) {
				size_t range = j->weights.size() * sizeof(decltype(j->weights)::value_type);
				flush_range(j->weights.data(), opts.stride, range);
				flush_range(j.get(), opts.stride, sizeof(JointAos));
			}

			size_t range = m.joints.size() * sizeof(decltype(m.joints)::value_type);
			flush_range(m.joints.data(), opts.stride, range);
		}
	}

	void run_animate_joints()
	{
		for (auto& m: models) {
			animate_joints(m.joints[0].get(),
						   m.input->joints_raw.data(),
						   m.input->frames[m.frame].values.data());
		}
	}

	void run_animate_weights()
	{
		for (auto& m: models) {
			animate_weights(m.joints[0].get());

			m.frame++;
			m.frame %= m.input->frames.size();
		}
	}

	void modify_structure(const BenchmarkOptions& opts)
	{
		size_t num_to_replace = opts.num_copies * TO_REMOVE_RATIO;
		for (size_t i = 0; i < num_to_replace; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				models.pop_front();
			}
		}

		for (size_t i = 0; i < num_to_replace; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				models.emplace_back(construct_model(INPUTS[j], counter));
			}
			counter++;
		}
	}
};

struct PooledJointsMixedWeights
{
	int frame;
	const ModelInput* input;
	std::vector<JointMixed> joints;
};

struct Bench_PooledJointsMixedWeights
{
	std::deque<PooledJointsMixedWeights> models;
	size_t counter = 0;

	PooledJointsMixedWeights construct_model(const ModelInput& input, size_t id)
	{
		PooledJointsMixedWeights model;
		model.input = &input;
		model.frame = id % input.frames.size();

		for (size_t i = 0; i < input.joints_raw.size(); i++) {
			JointMixed j;
			j.orient = input.joints_raw[i].base_quat;
			j.pos = input.joints_raw[i].base_pos;

			model.joints.emplace_back(std::move(j));
		}

		model.joints[0].parent = nullptr;
		model.joints[0].next = nullptr;
		for (size_t i = 1; i < model.joints.size(); i++) {
			model.joints[i].parent = &model.joints[input.joints_raw[i].parent];
			model.joints[i].next = nullptr;
			model.joints[i - 1].next = &model.joints[i];
		}

		for (size_t i = 0; i < input.weights_raw.size(); i++) {
			auto& joint = model.joints[input.weights_raw[i].joint];
			const auto& w = input.weights_raw[i];

			joint.weight_bias.push_back(w.bias);
			joint.weight_initial_pos.push_back(w.pos);

			joint.weight_pos.emplace_back();
		}

		return model;
	}

	void initialize_from_opts(BenchmarkOptions& opts)
	{
		for (size_t i = 0; i < opts.num_copies; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				models.emplace_back(construct_model(INPUTS[j], i));
			}
		}

		counter = opts.num_copies - 1;
	}

	void flush_from_cache(const BenchmarkOptions& opts)
	{
		for (const auto& m: models) {
			for (const auto& j: m.joints) {
				size_t range = j.weight_initial_pos.size() * sizeof(float);

				flush_range(j.weight_initial_pos.data(), opts.stride, 3 * range);
				flush_range(j.weight_pos.data(), opts.stride, 3 * range);

				flush_range(j.weight_bias.data(), opts.stride, range);
			}

			size_t range = m.joints.size() * sizeof(decltype(m.joints)::value_type);
			flush_range(m.joints.data(), opts.stride, range);
		}
	}

	void run_animate_joints()
	{
		for (auto& m: models) {
			animate_joints(&m.joints[0],
						   m.input->joints_raw.data(),
						   m.input->frames[m.frame].values.data());
		}
	}

	void run_animate_weights()
	{
		for (auto& m: models) {
			animate_weights(&m.joints[0]);

			m.frame++;
			m.frame %= m.input->frames.size();
		}
	}

	void modify_structure(const BenchmarkOptions& opts)
	{
		size_t num_to_replace = opts.num_copies * TO_REMOVE_RATIO;
		for (size_t i = 0; i < num_to_replace; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				models.pop_front();
			}
		}

		for (size_t i = 0; i < num_to_replace; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				models.emplace_back(construct_model(INPUTS[j], counter));
			}
			counter++;
		}
	}
};

struct ScatteredJointsMixedWeights
{
	int frame;
	const ModelInput* input;
	std::vector<std::unique_ptr<JointMixed>> joints;
};

struct Bench_ScatteredJointsMixedWeights
{
	std::deque<ScatteredJointsMixedWeights> models;
	size_t counter = 0;

	ScatteredJointsMixedWeights construct_model(const ModelInput& input, size_t id)
	{
		ScatteredJointsMixedWeights model;
		model.input = &input;
		model.frame = id % input.frames.size();

		for (size_t i = 0; i < input.joints_raw.size(); i++) {
			std::unique_ptr<JointMixed> j(new JointMixed);
			j->orient = input.joints_raw[i].base_quat;
			j->pos = input.joints_raw[i].base_pos;

			model.joints.emplace_back(std::move(j));
		}

		model.joints[0]->parent = nullptr;
		model.joints[0]->next = nullptr;
		for (size_t i = 1; i < model.joints.size(); i++) {
			model.joints[i]->parent = model.joints[input.joints_raw[i].parent].get();
			model.joints[i]->next = nullptr;
			model.joints[i - 1]->next = model.joints[i].get();
		}

		for (size_t i = 0; i < input.weights_raw.size(); i++) {
			auto& joint = model.joints[input.weights_raw[i].joint];
			const auto& w = input.weights_raw[i];

			joint->weight_bias.push_back(w.bias);
			joint->weight_initial_pos.push_back(w.pos);

			joint->weight_pos.emplace_back();
		}

		return model;
	}

	void initialize_from_opts(BenchmarkOptions& opts)
	{
		for (size_t i = 0; i < opts.num_copies; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				models.emplace_back(construct_model(INPUTS[j], i));
			}
		}

		counter = opts.num_copies - 1;
	}

	void flush_from_cache(const BenchmarkOptions& opts)
	{
		for (const auto& m: models) {
			for (const auto& j: m.joints) {
				size_t range = j->weight_initial_pos.size() * sizeof(float);

				flush_range(j->weight_initial_pos.data(), opts.stride, 3 * range);

				flush_range(j->weight_pos.data(), opts.stride, 3 * range);

				flush_range(j->weight_bias.data(), opts.stride, range);
			}

			size_t range = m.joints.size() * sizeof(decltype(m.joints)::value_type);
			flush_range(m.joints.data(), opts.stride, range);
		}
	}

	void run_animate_joints()
	{
		for (auto& m: models) {
			animate_joints(m.joints[0].get(),
						   m.input->joints_raw.data(),
						   m.input->frames[m.frame].values.data());
		}
	}

	void run_animate_weights()
	{
		for (auto& m: models) {
			animate_weights(m.joints[0].get());

			m.frame++;
			m.frame %= m.input->frames.size();
		}
	}

	void modify_structure(const BenchmarkOptions& opts)
	{
		size_t num_to_replace = opts.num_copies * TO_REMOVE_RATIO;
		for (size_t i = 0; i < num_to_replace; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				models.pop_front();
			}
		}

		for (size_t i = 0; i < num_to_replace; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				models.emplace_back(construct_model(INPUTS[j], counter));
			}
			counter++;
		}
	}
};

struct PooledJointsSoaWeights
{
	int frame;
	const ModelInput* input;
	std::vector<JointSoa> joints;
};

struct Bench_PooledJointsSoaWeights
{
	std::deque<PooledJointsSoaWeights> models;
	size_t counter = 0;

	PooledJointsSoaWeights construct_model(const ModelInput& input, size_t id)
	{
		PooledJointsSoaWeights model;
		model.input = &input;
		model.frame = id % input.frames.size();

		for (size_t i = 0; i < input.joints_raw.size(); i++) {
			JointSoa j;
			j.orient = input.joints_raw[i].base_quat;
			j.pos = input.joints_raw[i].base_pos;

			model.joints.emplace_back(std::move(j));
		}

		model.joints[0].parent = nullptr;
		model.joints[0].next = nullptr;
		for (size_t i = 1; i < model.joints.size(); i++) {
			model.joints[i].parent = &model.joints[input.joints_raw[i].parent];
			model.joints[i].next = nullptr;
			model.joints[i - 1].next = &model.joints[i];
		}

		for (size_t i = 0; i < input.weights_raw.size(); i++) {
			auto& joint = model.joints[input.weights_raw[i].joint];
			const auto& w = input.weights_raw[i];

			joint.weight_bias.push_back(w.bias);
			joint.weight_initial_posx.push_back(w.pos.x);
			joint.weight_initial_posy.push_back(w.pos.y);
			joint.weight_initial_posz.push_back(w.pos.z);

			joint.weight_posx.push_back(0);
			joint.weight_posy.push_back(0);
			joint.weight_posz.push_back(0);
		}

		return model;
	}

	void initialize_from_opts(BenchmarkOptions& opts)
	{
		for (size_t i = 0; i < opts.num_copies; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				models.emplace_back(construct_model(INPUTS[j], i));
			}
		}

		counter = opts.num_copies - 1;
	}

	void flush_from_cache(const BenchmarkOptions& opts)
	{
		for (const auto& m: models) {
			for (const auto& j: m.joints) {
				size_t range = j.weight_initial_posx.size() * sizeof(float);

				flush_range(j.weight_initial_posx.data(), opts.stride, range);
				flush_range(j.weight_initial_posy.data(), opts.stride, range);
				flush_range(j.weight_initial_posz.data(), opts.stride, range);

				flush_range(j.weight_posx.data(), opts.stride, range);
				flush_range(j.weight_posy.data(), opts.stride, range);
				flush_range(j.weight_posz.data(), opts.stride, range);

				flush_range(j.weight_bias.data(), opts.stride, range);
			}

			size_t range = m.joints.size() * sizeof(decltype(m.joints)::value_type);
			flush_range(m.joints.data(), opts.stride, range);
		}
	}

	void run_animate_joints()
	{
		for (auto& m: models) {
			animate_joints(&m.joints[0],
						   m.input->joints_raw.data(),
						   m.input->frames[m.frame].values.data());
		}
	}

	void run_animate_weights()
	{
		for (auto& m: models) {
			animate_weights(&m.joints[0]);

			m.frame++;
			m.frame %= m.input->frames.size();
		}
	}

	void modify_structure(const BenchmarkOptions& opts)
	{
		size_t num_to_replace = opts.num_copies * TO_REMOVE_RATIO;
		for (size_t i = 0; i < num_to_replace; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				models.pop_front();
			}
		}

		for (size_t i = 0; i < num_to_replace; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				models.emplace_back(construct_model(INPUTS[j], counter));
			}
			counter++;
		}
	}
};

struct ScatteredJointsSoaWeights
{
	int frame;
	const ModelInput* input;
	std::vector<std::unique_ptr<JointSoa>> joints;
};

struct Bench_ScatteredJointsSoaWeights
{
	std::deque<ScatteredJointsSoaWeights> models;
	size_t counter = 0;

	ScatteredJointsSoaWeights construct_model(const ModelInput& input, size_t id)
	{
		ScatteredJointsSoaWeights model;
		model.input = &input;
		model.frame = id % input.frames.size();

		for (size_t i = 0; i < input.joints_raw.size(); i++) {
			std::unique_ptr<JointSoa> j(new JointSoa);
			j->orient = input.joints_raw[i].base_quat;
			j->pos = input.joints_raw[i].base_pos;

			model.joints.emplace_back(std::move(j));
		}

		model.joints[0]->parent = nullptr;
		model.joints[0]->next = nullptr;
		for (size_t i = 1; i < model.joints.size(); i++) {
			model.joints[i]->parent = model.joints[input.joints_raw[i].parent].get();
			model.joints[i]->next = nullptr;
			model.joints[i - 1]->next = model.joints[i].get();
		}

		for (size_t i = 0; i < input.weights_raw.size(); i++) {
			auto& joint = model.joints[input.weights_raw[i].joint];
			const auto& w = input.weights_raw[i];

			joint->weight_bias.push_back(w.bias);
			joint->weight_initial_posx.push_back(w.pos.x);
			joint->weight_initial_posy.push_back(w.pos.y);
			joint->weight_initial_posz.push_back(w.pos.z);

			joint->weight_posx.push_back(0);
			joint->weight_posy.push_back(0);
			joint->weight_posz.push_back(0);
		}

		return model;
	}

	void initialize_from_opts(BenchmarkOptions& opts)
	{
		for (size_t i = 0; i < opts.num_copies; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				models.emplace_back(construct_model(INPUTS[j], i));
			}
		}

		counter = opts.num_copies - 1;
	}

	void flush_from_cache(const BenchmarkOptions& opts)
	{
		for (const auto& m: models) {
			for (const auto& j: m.joints) {
				size_t range = j->weight_initial_posx.size() * sizeof(float);

				flush_range(j->weight_initial_posx.data(), opts.stride, range);
				flush_range(j->weight_initial_posy.data(), opts.stride, range);
				flush_range(j->weight_initial_posz.data(), opts.stride, range);

				flush_range(j->weight_posx.data(), opts.stride, range);
				flush_range(j->weight_posy.data(), opts.stride, range);
				flush_range(j->weight_posz.data(), opts.stride, range);

				flush_range(j->weight_bias.data(), opts.stride, range);
			}

			size_t range = m.joints.size() * sizeof(decltype(m.joints)::value_type);
			flush_range(m.joints.data(), opts.stride, range);
		}
	}

	void run_animate_joints()
	{
		for (auto& m: models) {
			animate_joints(m.joints[0].get(),
						   m.input->joints_raw.data(),
						   m.input->frames[m.frame].values.data());
		}
	}

	void run_animate_weights()
	{
		for (auto& m: models) {
			animate_weights(m.joints[0].get());

			m.frame++;
			m.frame %= m.input->frames.size();
		}
	}

	void modify_structure(const BenchmarkOptions& opts)
	{
		size_t num_to_replace = opts.num_copies * TO_REMOVE_RATIO;
		for (size_t i = 0; i < num_to_replace; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				models.pop_front();
			}
		}

		for (size_t i = 0; i < num_to_replace; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				models.emplace_back(construct_model(INPUTS[j], counter));
			}
			counter++;
		}
	}
};

struct OnePoolJointsAosWeights
{
	int frame;
	const ModelInput* input;
	size_t idx_start;
	size_t idx_end;

};

struct Bench_OnePoolJointsAosWeights
{
	std::vector<OnePoolJointsAosWeights> models;
	std::vector<JointOnePool> joints;
	WeightPoolAos wp;

	size_t counter = 0;
	size_t joint_counter = 0;
	size_t weight_counter = 0;

	void construct_model(const ModelInput& input, size_t id)
	{
		models[counter].idx_start = joint_counter;
		models[counter].idx_end = joint_counter + input.joints_raw.size();
		models[counter].frame = id % input.frames.size();;
		models[counter].input = &input;
		for (size_t idx = 0; idx < input.joints_raw.size(); idx++) {
			size_t i = joint_counter + idx;
			joints[i].orient = input.joints_raw[idx].base_quat;
			joints[i].pos = input.joints_raw[idx].base_pos;
			joints[i].weights_start = weight_counter + input.joints_raw[idx].idx_start;
			joints[i].weights_end = weight_counter + input.joints_raw[idx].idx_end;
		}

		joints[joint_counter].parent = nullptr;
		joints[joint_counter].next = nullptr;
		for (size_t idx = 1; idx < input.joints_raw.size(); idx++) {
			size_t i = joint_counter + idx;
			joints[i].parent = &joints[joint_counter + input.joints_raw[idx].parent];
			joints[i].next = nullptr;
			joints[i - 1].next = &joints[i];
		}

		for (size_t idx = 0; idx < input.weights_raw.size(); idx++) {
			size_t i = weight_counter + idx;

			wp.weights[i].bias = input.weights_raw[idx].bias;
			wp.weights[i].initial_pos = input.weights_raw[idx].pos;
			wp.weights[i].pos.x = 0;
			wp.weights[i].pos.y = 0;
			wp.weights[i].pos.z = 0;
		}

		joint_counter += input.joints_raw.size();
		if (joint_counter == joints.size()) {
			joint_counter = 0;
		}

		weight_counter += input.weights_raw.size();
		if (weight_counter == wp.weights.size()) {
			weight_counter = 0;
		}

		counter++;
		if (counter == models.size()) {
			counter = 0;
		}
	}

	void initialize_from_opts(BenchmarkOptions& opts)
	{
		models.resize(opts.num_copies * NUM_MODELS);
		size_t num_joints = std::accumulate(INPUTS, INPUTS + NUM_MODELS, 0,
			[](size_t acc, const ModelInput& m) {
				return acc + m.joints_raw.size();
			});
		size_t num_weights = std::accumulate(INPUTS, INPUTS + NUM_MODELS, 0,
			[](size_t acc, const ModelInput& m) {
				return acc + m.weights_raw.size();
			});
		joints.resize(num_joints * models.size());
		wp.weights.resize(num_weights * models.size());

		for (size_t i = 0; i < opts.num_copies; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				construct_model(INPUTS[j], i);
			}
		}
	}

	void flush_from_cache(const BenchmarkOptions& opts)
	{
		size_t joints_range = joints.size() * sizeof(JointOnePool);
		flush_range(joints.data(), opts.stride, joints_range);

		size_t weights_range = wp.weights.size() * sizeof(WeightAos);
		flush_range(wp.weights.data(), opts.stride, weights_range);
	}

	void run_animate_joints()
	{
		for (auto& m: models) {
			JointOnePool* root = &joints[m.idx_start];
			animate_joints(root,
						   m.input->joints_raw.data(),
						   m.input->frames[m.frame].values.data());
		}
	}

	void run_animate_weights()
	{
		for (auto& m: models) {
			animate_weights(&joints[m.idx_start], wp);

			m.frame++;
			m.frame %= m.input->frames.size();
		}
	}

	void modify_structure(const BenchmarkOptions& opts)
	{
		size_t num_to_replace = opts.num_copies * TO_REMOVE_RATIO;
		for (size_t i = 0; i < num_to_replace; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				construct_model(INPUTS[j], i);
			}
		}
	}
};

struct OnePoolJointsMixedWeights
{
	int frame;
	const ModelInput* input;
	size_t idx_start;
	size_t idx_end;

};

struct Bench_OnePoolJointsMixedWeights
{
	std::vector<OnePoolJointsMixedWeights> models;
	std::vector<JointOnePool> joints;
	WeightPoolMixed wp;

	size_t counter = 0;
	size_t joint_counter = 0;
	size_t weight_counter = 0;

	void construct_model(const ModelInput& input, size_t id)
	{
		models[counter].idx_start = joint_counter;
		models[counter].idx_end = joint_counter + input.joints_raw.size();
		models[counter].frame = id % input.frames.size();;
		models[counter].input = &input;
		for (size_t idx = 0; idx < input.joints_raw.size(); idx++) {
			size_t i = joint_counter + idx;
			joints[i].orient = input.joints_raw[idx].base_quat;
			joints[i].pos = input.joints_raw[idx].base_pos;
			joints[i].weights_start = weight_counter + input.joints_raw[idx].idx_start;
			joints[i].weights_end = weight_counter + input.joints_raw[idx].idx_end;
		}

		joints[joint_counter].parent = nullptr;
		joints[joint_counter].next = nullptr;
		for (size_t idx = 1; idx < input.joints_raw.size(); idx++) {
			size_t i = joint_counter + idx;
			joints[i].parent = &joints[joint_counter + input.joints_raw[idx].parent];
			joints[i].next = nullptr;
			joints[i - 1].next = &joints[i];
		}

		for (size_t idx = 0; idx < input.weights_raw.size(); idx++) {
			size_t i = weight_counter + idx;

			wp.bias[i] = input.weights_raw[idx].bias;
			wp.initial_pos[i] = input.weights_raw[idx].pos;
			wp.pos[i].x = 0;
			wp.pos[i].y = 0;
			wp.pos[i].z = 0;
		}

		joint_counter += input.joints_raw.size();
		if (joint_counter == joints.size()) {
			joint_counter = 0;
		}

		weight_counter += input.weights_raw.size();
		if (weight_counter == wp.pos.size()) {
			weight_counter = 0;
		}

		counter++;
		if (counter == models.size()) {
			counter = 0;
		}
	}

	void initialize_from_opts(BenchmarkOptions& opts)
	{
		models.resize(opts.num_copies * NUM_MODELS);
		size_t num_joints = std::accumulate(INPUTS, INPUTS + NUM_MODELS, 0,
			[](size_t acc, const ModelInput& m) {
				return acc + m.joints_raw.size();
			});
		size_t num_weights = std::accumulate(INPUTS, INPUTS + NUM_MODELS, 0,
			[](size_t acc, const ModelInput& m) {
				return acc + m.weights_raw.size();
			});
		joints.resize(num_joints * models.size());

		wp.pos.resize(num_weights * models.size());
		wp.bias.resize(num_weights * models.size());
		wp.initial_pos.resize(num_weights * models.size());

		for (size_t i = 0; i < opts.num_copies; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				construct_model(INPUTS[j], i);
			}
		}
	}

	void flush_from_cache(const BenchmarkOptions& opts)
	{
		size_t joints_range = joints.size() * sizeof(JointOnePool);
		flush_range(joints.data(), opts.stride, joints_range);

		size_t weights_range = wp.pos.size() * sizeof(float);
		flush_range(wp.pos.data(), opts.stride, 3 * weights_range);
		flush_range(wp.initial_pos.data(), opts.stride, 3 * weights_range);
		flush_range(wp.bias.data(), opts.stride, weights_range);
	}

	void run_animate_joints()
	{
		for (auto& m: models) {
			JointOnePool* root = &joints[m.idx_start];
			animate_joints(root,
						   m.input->joints_raw.data(),
						   m.input->frames[m.frame].values.data());
		}
	}

	void run_animate_weights()
	{
		for (auto& m: models) {
			animate_weights(&joints[m.idx_start], wp);

			m.frame++;
			m.frame %= m.input->frames.size();
		}
	}

	void modify_structure(const BenchmarkOptions& opts)
	{
		size_t num_to_replace = opts.num_copies * TO_REMOVE_RATIO;
		for (size_t i = 0; i < num_to_replace; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				construct_model(INPUTS[j], i);
			}
		}
	}
};

struct OnePoolJointsSoaWeights
{
	int frame;
	const ModelInput* input;
	size_t idx_start;
	size_t idx_end;

};

struct Bench_OnePoolJointsSoaWeights
{
	std::vector<OnePoolJointsSoaWeights> models;
	std::vector<JointOnePool> joints;
	WeightPoolSoa wp;

	size_t counter = 0;
	size_t joint_counter = 0;
	size_t weight_counter = 0;

	void construct_model(const ModelInput& input, size_t id)
	{
		models[counter].idx_start = joint_counter;
		models[counter].idx_end = joint_counter + input.joints_raw.size();
		models[counter].frame = id % input.frames.size();;
		models[counter].input = &input;
		for (size_t idx = 0; idx < input.joints_raw.size(); idx++) {
			size_t i = joint_counter + idx;
			joints[i].orient = input.joints_raw[idx].base_quat;
			joints[i].pos = input.joints_raw[idx].base_pos;
			joints[i].weights_start = weight_counter + input.joints_raw[idx].idx_start;
			joints[i].weights_end = weight_counter + input.joints_raw[idx].idx_end;
		}

		joints[joint_counter].parent = nullptr;
		joints[joint_counter].next = nullptr;
		for (size_t idx = 1; idx < input.joints_raw.size(); idx++) {
			size_t i = joint_counter + idx;
			joints[i].parent = &joints[joint_counter + input.joints_raw[idx].parent];
			joints[i].next = nullptr;
			joints[i - 1].next = &joints[i];
		}

		for (size_t idx = 0; idx < input.weights_raw.size(); idx++) {
			size_t i = weight_counter + idx;

			wp.bias[i] = input.weights_raw[idx].bias;
			wp.initial_posx[i] = input.weights_raw[idx].pos.x;
			wp.initial_posy[i] = input.weights_raw[idx].pos.y;
			wp.initial_posz[i] = input.weights_raw[idx].pos.z;

			wp.posx[i] = 0;
			wp.posy[i] = 0;
			wp.posz[i] = 0;
		}

		joint_counter += input.joints_raw.size();
		if (joint_counter == joints.size()) {
			joint_counter = 0;
		}

		weight_counter += input.weights_raw.size();
		if (weight_counter == wp.posx.size()) {
			weight_counter = 0;
		}

		counter++;
		if (counter == models.size()) {
			counter = 0;
		}
	}

	void initialize_from_opts(BenchmarkOptions& opts)
	{
		models.resize(opts.num_copies * NUM_MODELS);
		size_t num_joints = std::accumulate(INPUTS, INPUTS + NUM_MODELS, 0,
			[](size_t acc, const ModelInput& m) {
				return acc + m.joints_raw.size();
			});
		size_t num_weights = std::accumulate(INPUTS, INPUTS + NUM_MODELS, 0,
			[](size_t acc, const ModelInput& m) {
				return acc + m.weights_raw.size();
			});
		joints.resize(num_joints * models.size());

		wp.posx.resize(num_weights * models.size());
		wp.posy.resize(num_weights * models.size());
		wp.posz.resize(num_weights * models.size());

		wp.initial_posx.resize(num_weights * models.size());
		wp.initial_posy.resize(num_weights * models.size());
		wp.initial_posz.resize(num_weights * models.size());

		wp.bias.resize(num_weights * models.size());

		for (size_t i = 0; i < opts.num_copies; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				construct_model(INPUTS[j], i);
			}
		}
	}

	void flush_from_cache(const BenchmarkOptions& opts)
	{
		size_t joints_range = joints.size() * sizeof(JointOnePool);
		flush_range(joints.data(), opts.stride, joints_range);

		size_t weights_range = wp.posx.size() * sizeof(float);
		flush_range(wp.posx.data(), opts.stride, weights_range);
		flush_range(wp.posy.data(), opts.stride, weights_range);
		flush_range(wp.posz.data(), opts.stride, weights_range);

		flush_range(wp.initial_posx.data(), opts.stride, weights_range);
		flush_range(wp.initial_posy.data(), opts.stride, weights_range);
		flush_range(wp.initial_posz.data(), opts.stride, weights_range);

		flush_range(wp.bias.data(), opts.stride, weights_range);
	}

	void run_animate_joints()
	{
		for (auto& m: models) {
			JointOnePool* root = &joints[m.idx_start];
			animate_joints(root,
						   m.input->joints_raw.data(),
						   m.input->frames[m.frame].values.data());
		}
	}

	void run_animate_weights()
	{
		for (auto& m: models) {
			animate_weights(&joints[m.idx_start], wp);

			m.frame++;
			m.frame %= m.input->frames.size();
		}
	}

	void modify_structure(const BenchmarkOptions& opts)
	{
		size_t num_to_replace = opts.num_copies * TO_REMOVE_RATIO;
		for (size_t i = 0; i < num_to_replace; i++) {
			for (size_t j = 0; j < NUM_MODELS; j++) {
				construct_model(INPUTS[j], i);
			}
		}
	}
};

template<typename DataStructure>
void BM_Template(benchmark::State& state)
{
	if (!initialized) {
		ModelInput input0 = parse_files("../resources/hellknight.md5mesh", "../resources/hellknight.md5anim");
		ModelInput input1 = parse_files("../resources/boblampclean.md5mesh", "../resources/boblampclean.md5anim");
		INPUTS[0] = std::move(input0);
		INPUTS[1] = std::move(input1);
		initialized = true;
	}

	auto opts = from_args(state);
	DataStructure structure;
	structure.initialize_from_opts(opts);

	for (auto _: state) {
		state.PauseTiming();
		if (opts.flush_range) {
			structure.flush_from_cache(opts);
		}
		structure.run_animate_joints();

		state.ResumeTiming();
		benchmark::ClobberMemory();
		structure.run_animate_weights();
		state.PauseTiming();

		benchmark::ClobberMemory();
		structure.modify_structure(opts);
		state.ResumeTiming();
	}

	size_t num_weights = 0;
	for (size_t i = 0; i < NUM_MODELS; i++) {
		num_weights += INPUTS[i].weights_raw.size();
	}

	state.SetComplexityN(opts.num_copies);
	auto items = state.iterations() * opts.num_copies * num_weights;
	state.SetItemsProcessed(items);
	state.SetBytesProcessed(6 * sizeof(float) * items);
}

double max_of_vector(const std::vector<double>& v)
{
	return *(std::max_element(v.begin(), v.end()));
}

double min_of_vector(const std::vector<double>& v)
{
	return *(std::min_element(v.begin(), v.end()));
}

BENCHMARK_TEMPLATE(BM_Template, Bench_ScatteredJointsAosWeights)
	->ArgNames({"NumDuplicates", "FlushCache"})
	->Apply(CustomArgs)
	->Complexity(benchmark::oN)
	->ComputeStatistics("min", min_of_vector)
	->ComputeStatistics("max", max_of_vector);

BENCHMARK_TEMPLATE(BM_Template, Bench_PooledJointsAosWeights)
	->ArgNames({"NumDuplicates", "FlushCache"})
	->Apply(CustomArgs)
	->Complexity(benchmark::oN)
	->ComputeStatistics("min", min_of_vector)
	->ComputeStatistics("max", max_of_vector);

BENCHMARK_TEMPLATE(BM_Template, Bench_ScatteredJointsMixedWeights)
	->ArgNames({"NumDuplicates", "FlushCache"})
	->Apply(CustomArgs)
	->Complexity(benchmark::oN)
	->ComputeStatistics("min", min_of_vector)
	->ComputeStatistics("max", max_of_vector);

BENCHMARK_TEMPLATE(BM_Template, Bench_PooledJointsMixedWeights)
	->ArgNames({"NumDuplicates", "FlushCache"})
	->Apply(CustomArgs)
	->Complexity(benchmark::oN)
	->ComputeStatistics("min", min_of_vector)
	->ComputeStatistics("max", max_of_vector);

BENCHMARK_TEMPLATE(BM_Template, Bench_ScatteredJointsSoaWeights)
	->ArgNames({"NumDuplicates", "FlushCache"})
	->Apply(CustomArgs)
	->Complexity(benchmark::oN)
	->ComputeStatistics("min", min_of_vector)
	->ComputeStatistics("max", max_of_vector);

BENCHMARK_TEMPLATE(BM_Template, Bench_PooledJointsSoaWeights)
	->ArgNames({"NumDuplicates", "FlushCache"})
	->Apply(CustomArgs)
	->Complexity(benchmark::oN)
	->ComputeStatistics("min", min_of_vector)
	->ComputeStatistics("max", max_of_vector);

BENCHMARK_TEMPLATE(BM_Template, Bench_OnePoolJointsAosWeights)
	->ArgNames({"NumDuplicates", "FlushCache"})
	->Apply(CustomArgs)
	->Complexity(benchmark::oN)
	->ComputeStatistics("min", min_of_vector)
	->ComputeStatistics("max", max_of_vector);

BENCHMARK_TEMPLATE(BM_Template, Bench_OnePoolJointsMixedWeights)
	->ArgNames({"NumDuplicates", "FlushCache"})
	->Apply(CustomArgs)
	->Complexity(benchmark::oN)
	->ComputeStatistics("min", min_of_vector)
	->ComputeStatistics("max", max_of_vector);

BENCHMARK_TEMPLATE(BM_Template, Bench_OnePoolJointsSoaWeights)
	->ArgNames({"NumDuplicates", "FlushCache"})
	->Apply(CustomArgs)
	->Complexity(benchmark::oN)
	->ComputeStatistics("min", min_of_vector)
	->ComputeStatistics("max", max_of_vector);

BENCHMARK_MAIN();
