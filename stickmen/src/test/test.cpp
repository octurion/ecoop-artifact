#include "anim/anim.h"

#include "test/dataset.h"

#include <benchmark/benchmark.h>

#include <unordered_map>
#include <utility>

#include <stdio.h>
#include <unistd.h>

enum class DatasetType {
	Dataset1,
	Dataset2,
	DatasetCount
};

struct BenchmarkOptions
{
	const JointData* joint_data;
	size_t joint_count;

	const WeightData* weight_data;
	size_t weight_count;

	size_t weight_replication_count;
	size_t stride;
	bool flush_range;
};

static BenchmarkOptions from_args(benchmark::State& state)
{
	BenchmarkOptions opt;

	opt.weight_replication_count = state.range(0);

	switch (static_cast<DatasetType>(state.range(1))) {
	case DatasetType::Dataset1:
		opt.joint_data = JOINTS1;
		opt.joint_count = NUM_JOINTS1;

		opt.weight_data = WEIGHTS1;
		opt.weight_count = NUM_WEIGHTS1;
		break;

	case DatasetType::Dataset2:
		opt.joint_data = JOINTS2;
		opt.joint_count = NUM_JOINTS2;

		opt.weight_data = WEIGHTS2;
		opt.weight_count = NUM_WEIGHTS2;
		break;

	default:
		abort();
	}

	opt.stride = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);

	return opt;
}

static void CustomArgs(benchmark::internal::Benchmark* b)
{
	for (int i = 1; i <= 10; i++) {
		b->Args({i * 32, (int)DatasetType::Dataset1});
	}
}

static void flush_range(const void* ptr, size_t stride, size_t len)
{
	size_t flush_count = (len + stride) & (stride - 1);

	for (size_t i = 0; i < flush_count; i += stride) {
		_mm_clflush((const char*) ptr + i * stride);
	}
}

struct ScatteredOwningAos
{
	JointScatteredOwning root;
	std::vector<std::pair<WeightParentless*, vec3>> initial;

	void initialize_from_opts(BenchmarkOptions& opts)
	{
		std::unordered_map<int, JointScatteredOwning*> joint_map;
		joint_map[0] = &root;

		root.orient_relative = quat::from_xyz(
			opts.joint_data[0].qx, opts.joint_data[0].qy, opts.joint_data[0].qz);

		root.orient_absolute.w = 1;
		root.orient_absolute.x = 0;
		root.orient_absolute.y = 0;
		root.orient_absolute.z = 0;

		for (size_t i = 1; i < opts.joint_count; i++) {
			auto& parent = joint_map[opts.joint_data[i].parent];

			std::unique_ptr<JointScatteredOwning> new_joint(new JointScatteredOwning);
			new_joint->orient_relative = quat::from_xyz(
				opts.joint_data[i].qx, opts.joint_data[i].qy, opts.joint_data[i].qz);

			new_joint->orient_absolute.w = 1;
			new_joint->orient_absolute.x = 0;
			new_joint->orient_absolute.y = 0;
			new_joint->orient_absolute.z = 0;

			joint_map[i] = new_joint.get();
			parent->children.emplace_back(std::move(new_joint));
		}

		for (size_t c = 0; c < opts.weight_replication_count; c++) {
			for (size_t i = 0; i < opts.weight_count; i++) {
				auto& joint = joint_map[opts.weight_data[i].joint];

				WeightParentless new_weight;

				new_weight.w = opts.weight_data[i].w;

				new_weight.pos.x = opts.weight_data[i].x;
				new_weight.pos.y = opts.weight_data[i].y;
				new_weight.pos.z = opts.weight_data[i].z;

				joint->weights.push_back(new_weight);
			}
		}

		for (auto& j: joint_map) {
			for (auto& e: j.second->weights) {
				initial.emplace_back(&e, e.pos);
			}
		}
	}

	void flush_from_cache(const BenchmarkOptions& opts)
	{
		flush_from_cache(root, opts.stride);
	}

	static void flush_from_cache(JointScatteredOwning& node, size_t stride)
	{
		for (auto& e: node.children) {
			flush_from_cache(*e, stride);
		}

		flush_range(node.weights.data(),
					stride,
					node.weights.size() * sizeof(decltype(root.weights)::value_type));

		flush_range(node.children.data(),
					stride,
					node.children.size() * sizeof(decltype(root.children)::value_type));

		flush_range(&node, stride, sizeof(root));
	}

	void animate_joints()
	{
		animate_joints_scattered_owning(root);
	}

	void animate_weights()
	{
		animate_weights_scattered_parentless(root);
	}

	void reset_weights()
	{
		for (auto& e: initial) {
			e.first->pos = e.second;
		}
	}
};

struct PooledOwningAos
{
	std::vector<JointPooledOwning> joints;
	std::vector<std::pair<WeightParentless*, vec3>> initial;

	void initialize_from_opts(BenchmarkOptions& opts)
	{
		joints.resize(opts.joint_count);
		auto& root = joints[0];
		root.orient_relative = quat::from_xyz(
			opts.joint_data[0].qx, opts.joint_data[0].qy, opts.joint_data[0].qz);

		root.orient_absolute.w = 1;
		root.orient_absolute.x = 0;
		root.orient_absolute.y = 0;
		root.orient_absolute.z = 0;

		for (size_t i = 1; i < opts.joint_count; i++) {
			auto& new_joint = joints[i];
			new_joint.orient_relative = quat::from_xyz(
				opts.joint_data[i].qx, opts.joint_data[i].qy, opts.joint_data[i].qz);

			new_joint.orient_absolute.w = 1;
			new_joint.orient_absolute.x = 0;
			new_joint.orient_absolute.y = 0;
			new_joint.orient_absolute.z = 0;

			new_joint.parent = &joints[opts.joint_data[i].parent];
		}

		for (size_t c = 0; c < opts.weight_replication_count; c++) {
			for (size_t i = 0; i < opts.weight_count; i++) {
				auto& joint = joints[opts.weight_data[i].joint];
				WeightParentless new_weight;

				new_weight.w = opts.weight_data[i].w;
				new_weight.pos.x = opts.weight_data[i].x;
				new_weight.pos.y = opts.weight_data[i].y;
				new_weight.pos.z = opts.weight_data[i].z;

				joint.weights.push_back(new_weight);
			}
		}

		for (auto& j: joints) {
			for (auto& e: j.weights) {
				initial.emplace_back(&e, e.pos);
			}
		}
	}

	void flush_from_cache(const BenchmarkOptions& opts)
	{
		for (auto& e: joints) {
			flush_range(e.weights.data(), opts.stride,
						e.weights.size() * sizeof(decltype(e.weights)::value_type));
		}
		flush_range(joints.data(), opts.stride,
					joints.size() * sizeof(decltype(joints)::value_type));
	}

	void animate_joints()
	{
		animate_joints_pooled_owning(joints);
	}

	void animate_weights()
	{
		animate_weights_pooled_parentless(joints);
	}

	void reset_weights()
	{
		for (auto& e: initial) {
			e.first->pos = e.second;
		}
	}
};

struct ScatteredOwnerlessAos
{
	JointScatteredOwnerless root;
	std::vector<WeightScatteredParent> weights;
	std::vector<vec3> initial;

	void initialize_from_opts(BenchmarkOptions& opts)
	{
		std::unordered_map<int, JointScatteredOwnerless*> joint_map;

		joint_map[0] = &root;

		weights.reserve(opts.weight_count);

		root.orient_relative = quat::from_xyz(
			opts.joint_data[0].qx, opts.joint_data[0].qy, opts.joint_data[0].qz);

		root.orient_absolute.w = 1;
		root.orient_absolute.x = 0;
		root.orient_absolute.y = 0;
		root.orient_absolute.z = 0;

		for (size_t i = 1; i < opts.joint_count; i++) {
			auto& parent = joint_map[opts.joint_data[i].parent];

			std::unique_ptr<JointScatteredOwnerless> new_joint(new JointScatteredOwnerless);
			new_joint->orient_relative = quat::from_xyz(
				opts.joint_data[i].qx, opts.joint_data[i].qy, opts.joint_data[i].qz);

			new_joint->orient_absolute.w = 1;
			new_joint->orient_absolute.x = 0;
			new_joint->orient_absolute.y = 0;
			new_joint->orient_absolute.z = 0;

			joint_map[i] = new_joint.get();
			parent->children.emplace_back(std::move(new_joint));
		}

		for (size_t c = 0; c < opts.weight_replication_count; c++) {
			for (size_t i = 0; i < opts.weight_count; i++) {
				auto* joint = joint_map[opts.weight_data[i].joint];

				WeightScatteredParent new_weight;

				new_weight.w = opts.weight_data[i].w;

				new_weight.pos.x = opts.weight_data[i].x;
				new_weight.pos.y = opts.weight_data[i].y;
				new_weight.pos.z = opts.weight_data[i].z;
				new_weight.joint = joint;

				weights.push_back(new_weight);
				initial.push_back(new_weight.pos);
			}
		}
	}

	void flush_from_cache(const BenchmarkOptions& opts)
	{
		flush_from_cache(root, opts.stride);
		flush_range(weights.data(), opts.stride,
					weights.size() * sizeof(decltype(weights)::value_type));
	}

	static void flush_from_cache(JointScatteredOwnerless& root, size_t stride)
	{
		for (auto& e: root.children) {
			flush_from_cache(*e, stride);
		}

		flush_range(root.children.data(),
					stride,
					root.children.size() * sizeof(decltype(root.children)::value_type));

		flush_range(&root, stride, sizeof(root));
	}

	void animate_joints()
	{
		animate_joints_scattered_ownerless(root);
	}

	void animate_weights()
	{
		animate_weights_scattered_parent(weights);
	}

	void reset_weights()
	{
		for (size_t i = 0; i < initial.size(); i++) {
			weights[i].pos = initial[i];
		}
	}
};

struct PooledOwnerlessAos
{
	std::vector<JointPooledOwnerless> joints;
	std::vector<WeightPooledParent> weights;
	std::vector<vec3> initial;

	void initialize_from_opts(BenchmarkOptions& opts)
	{
		joints.resize(opts.joint_count);
		weights.reserve(opts.weight_count);

		auto& root = joints[0];
		root.orient_relative = quat::from_xyz(
			opts.joint_data[0].qx, opts.joint_data[0].qy, opts.joint_data[0].qz);

		root.orient_absolute.w = 1;
		root.orient_absolute.x = 0;
		root.orient_absolute.y = 0;
		root.orient_absolute.z = 0;

		for (size_t i = 1; i < opts.joint_count; i++) {
			auto& new_joint = joints[i];
			new_joint.orient_relative = quat::from_xyz(
				opts.joint_data[i].qx, opts.joint_data[i].qy, opts.joint_data[i].qz);

			new_joint.orient_absolute.w = 1;
			new_joint.orient_absolute.x = 0;
			new_joint.orient_absolute.y = 0;
			new_joint.orient_absolute.z = 0;

			new_joint.parent = &joints[opts.joint_data[i].parent];
		}

		for (size_t c = 0; c < opts.weight_replication_count; c++) {
			for (size_t i = 0; i < opts.weight_count; i++) {
				auto* joint = &joints[opts.weight_data[i].joint];
				WeightPooledParent new_weight;

				new_weight.w = opts.weight_data[i].w;
				new_weight.pos.x = opts.weight_data[i].x;
				new_weight.pos.y = opts.weight_data[i].y;
				new_weight.pos.z = opts.weight_data[i].z;
				new_weight.joint = joint;

				weights.push_back(new_weight);
				initial.push_back(new_weight.pos);
			}
		}
	}

	void flush_from_cache(const BenchmarkOptions& opts)
	{
		flush_range(joints.data(),
					opts.stride,
					joints.size() * sizeof(decltype(joints)::value_type));
		flush_range(weights.data(),
					opts.stride,
					weights.size() * sizeof(decltype(weights)::value_type));
	}

	void animate_joints()
	{
		animate_joints_pooled_ownerless(joints);
	}

	void animate_weights()
	{
		animate_weights_pooled_parent(weights);
	}

	void reset_weights()
	{
		for (size_t i = 0; i < initial.size(); i++) {
			weights[i].pos = initial[i];
		}
	}
};

struct ScatteredOwningSoa
{
	struct ResetInfo {
		JointScatteredOwningSoa* joint;
		size_t idx;
		vec3 pos;

		ResetInfo(JointScatteredOwningSoa* joint, size_t idx, vec3 pos)
			: joint(joint)
			, idx(idx)
			, pos(pos)
		{
		}
	};

	JointScatteredOwningSoa root;

	std::vector<ResetInfo> initial;

	void initialize_from_opts(BenchmarkOptions& opts)
	{
		std::unordered_map<int, JointScatteredOwningSoa*> joint_map;
		joint_map[0] = &root;

		root.orient_relative = quat::from_xyz(
			opts.joint_data[0].qx, opts.joint_data[0].qy, opts.joint_data[0].qz);

		root.orient_absolute.w = 1;
		root.orient_absolute.x = 0;
		root.orient_absolute.y = 0;
		root.orient_absolute.z = 0;

		for (size_t i = 1; i < opts.joint_count; i++) {
			auto& parent = joint_map[opts.joint_data[i].parent];

			std::unique_ptr<JointScatteredOwningSoa> new_joint(new JointScatteredOwningSoa);
			new_joint->orient_relative = quat::from_xyz(
				opts.joint_data[i].qx, opts.joint_data[i].qy, opts.joint_data[i].qz);

			new_joint->orient_absolute.w = 1;
			new_joint->orient_absolute.x = 0;
			new_joint->orient_absolute.y = 0;
			new_joint->orient_absolute.z = 0;

			joint_map[i] = new_joint.get();
			parent->children.emplace_back(std::move(new_joint));
		}

		for (size_t c = 0; c < opts.weight_replication_count; c++) {
			for (size_t i = 0; i < opts.weight_count; i++) {
				auto& joint = joint_map[opts.weight_data[i].joint];

				joint->weights.w.push_back(opts.weight_data[i].w);
				joint->weights.posx.push_back(opts.weight_data[i].x);
				joint->weights.posy.push_back(opts.weight_data[i].y);
				joint->weights.posz.push_back(opts.weight_data[i].z);

				vec3 v;
				v.x = opts.weight_data[i].x;
				v.y = opts.weight_data[i].y;
				v.z = opts.weight_data[i].z;

				initial.emplace_back(joint, joint->weights.w.size() - 1, v);
			}
		}
	}

	static void flush_from_cache(JointScatteredOwningSoa& root, size_t stride)
	{
		for (auto& e: root.children) {
			flush_from_cache(*e, stride);
		}

		auto& weights = root.weights;

		flush_range(weights.w.data(),
					stride,
					weights.w.size() * sizeof(decltype(weights.w)::value_type));
		flush_range(weights.posx.data(),
					stride,
					weights.posx.size() * sizeof(decltype(weights.posx)::value_type));
		flush_range(weights.posy.data(),
					stride,
					weights.posy.size() * sizeof(decltype(weights.posy)::value_type));
		flush_range(weights.posz.data(),
					stride,
					weights.posz.size() * sizeof(decltype(weights.posz)::value_type));

		flush_range(root.children.data(),
					stride,
					root.children.size() * sizeof(decltype(root.children)::value_type));

		flush_range(&root, stride, sizeof(root));
	}

	void flush_from_cache(const BenchmarkOptions& opts)
	{
		flush_from_cache(root, opts.stride);
	}

	void animate_joints()
	{
		animate_joints_scattered_owning_soa(root);
	}

	void animate_weights()
	{
		animate_weights_scattered_parentless_soa(root);
	}

	void reset_weights()
	{
		for (auto& e: initial) {
			e.joint->weights.posx[e.idx] = e.pos.x;
			e.joint->weights.posy[e.idx] = e.pos.y;
			e.joint->weights.posz[e.idx] = e.pos.z;
		}
	}
};

struct PooledOwningSoa
{
	struct ResetInfo {
		JointPooledOwningSoa* joint;
		size_t idx;
		vec3 pos;

		ResetInfo(JointPooledOwningSoa* joint, size_t idx, vec3 pos)
			: joint(joint)
			, idx(idx)
			, pos(pos)
		{
		}
	};

	std::vector<JointPooledOwningSoa> joints;

	std::vector<ResetInfo> initial;

	void initialize_from_opts(BenchmarkOptions& opts)
	{
		joints.resize(opts.joint_count);
		auto& root = joints[0];
		root.orient_relative = quat::from_xyz(
			opts.joint_data[0].qx, opts.joint_data[0].qy, opts.joint_data[0].qz);

		root.orient_absolute.w = 1;
		root.orient_absolute.x = 0;
		root.orient_absolute.y = 0;
		root.orient_absolute.z = 0;

		for (size_t i = 1; i < opts.joint_count; i++) {
			auto& new_joint = joints[i];
			new_joint.orient_relative = quat::from_xyz(
				opts.joint_data[i].qx, opts.joint_data[i].qy, opts.joint_data[i].qz);

			new_joint.orient_absolute.w = 1;
			new_joint.orient_absolute.x = 0;
			new_joint.orient_absolute.y = 0;
			new_joint.orient_absolute.z = 0;

			new_joint.parent = &joints[opts.joint_data[i].parent];
		}

		for (size_t c = 0; c < opts.weight_replication_count; c++) {
			for (size_t i = 0; i < opts.weight_count; i++) {
				auto& joint = joints[opts.weight_data[i].joint];

				joint.weights.w.push_back(opts.weight_data[i].w);

				joint.weights.posx.push_back(opts.weight_data[i].x);
				joint.weights.posy.push_back(opts.weight_data[i].y);
				joint.weights.posz.push_back(opts.weight_data[i].z);

				vec3 v;
				v.x = opts.weight_data[i].x;
				v.y = opts.weight_data[i].y;
				v.z = opts.weight_data[i].z;

				initial.emplace_back(&joint, joint.weights.w.size() - 1, v);
			}
		}
	}

	void flush_from_cache(const BenchmarkOptions& opts)
	{
		for (auto& e: joints) {
			flush_range(e.weights.w.data(), opts.stride,
						e.weights.w.size() * sizeof(decltype(e.weights.w)::value_type));

			flush_range(e.weights.posx.data(), opts.stride,
						e.weights.posx.size() * sizeof(decltype(e.weights.posx)::value_type));
			flush_range(e.weights.posy.data(), opts.stride,
						e.weights.posy.size() * sizeof(decltype(e.weights.posy)::value_type));
			flush_range(e.weights.posz.data(), opts.stride,
						e.weights.posz.size() * sizeof(decltype(e.weights.posz)::value_type));
		}
		flush_range(joints.data(), opts.stride,
					joints.size() * sizeof(decltype(joints)::value_type));
	}

	void animate_joints()
	{
		animate_joints_pooled_owning_soa(joints);
	}

	void animate_weights()
	{
		animate_weights_pooled_parentless_soa(joints);
	}

	void reset_weights()
	{
		for (auto& e: initial) {
			e.joint->weights.posx[e.idx] = e.pos.x;
			e.joint->weights.posy[e.idx] = e.pos.y;
			e.joint->weights.posz[e.idx] = e.pos.z;
		}
	}
};

struct ScatteredOwnerlessSoa
{
	JointScatteredOwnerlessSoa root;
	WeightScatteredParentSoa weights;
	std::vector<vec3> initial;

	void initialize_from_opts(BenchmarkOptions& opts)
	{
		std::unordered_map<int, JointScatteredOwnerlessSoa*> joint_map;

		joint_map[0] = &root;

		weights.joint.reserve(opts.weight_count);
		weights.w.reserve(opts.weight_count);
		weights.posx.reserve(opts.weight_count);
		weights.posy.reserve(opts.weight_count);
		weights.posz.reserve(opts.weight_count);

		root.orient_relative = quat::from_xyz(
			opts.joint_data[0].qx, opts.joint_data[0].qy, opts.joint_data[0].qz);

		root.orient_absolute.w = 1;
		root.orient_absolute.x = 0;
		root.orient_absolute.y = 0;
		root.orient_absolute.z = 0;

		for (size_t i = 1; i < opts.joint_count; i++) {
			auto& parent = joint_map[opts.joint_data[i].parent];

			std::unique_ptr<JointScatteredOwnerlessSoa> new_joint(new JointScatteredOwnerlessSoa);
			new_joint->orient_relative = quat::from_xyz(
				opts.joint_data[i].qx, opts.joint_data[i].qy, opts.joint_data[i].qz);

			new_joint->orient_absolute.w = 1;
			new_joint->orient_absolute.x = 0;
			new_joint->orient_absolute.y = 0;
			new_joint->orient_absolute.z = 0;

			joint_map[i] = new_joint.get();
			parent->children.emplace_back(std::move(new_joint));
		}

		for (size_t c = 0; c < opts.weight_replication_count; c++) {
			for (size_t i = 0; i < opts.weight_count; i++) {
				auto* joint = joint_map[opts.weight_data[i].joint];

				weights.joint.push_back(joint);
				weights.w.push_back(opts.weight_data[i].w);

				weights.posx.push_back(opts.weight_data[i].x);
				weights.posy.push_back(opts.weight_data[i].y);
				weights.posz.push_back(opts.weight_data[i].z);

				vec3 v;
				v.x = opts.weight_data[i].x;
				v.y = opts.weight_data[i].y;
				v.z = opts.weight_data[i].z;

				initial.push_back(v);
			}
		}
	}

	static void flush_from_cache(JointScatteredOwnerlessSoa& root, size_t stride)
	{
		for (auto& e: root.children) {
			flush_from_cache(*e, stride);
		}

		flush_range(root.children.data(),
					stride,
					root.children.size() * sizeof(decltype(root.children)::value_type));

		flush_range(&root, stride, sizeof(root));
	}

	void flush_from_cache(const BenchmarkOptions& opts)
	{
			flush_from_cache(root, opts.stride);
			flush_range(weights.joint.data(), opts.stride,
						weights.joint.size() * sizeof(decltype(weights.joint)::value_type));
			flush_range(weights.w.data(), opts.stride,
						weights.w.size() * sizeof(decltype(weights.w)::value_type));
			flush_range(weights.posx.data(), opts.stride,
						weights.posx.size() * sizeof(decltype(weights.posx)::value_type));
			flush_range(weights.posy.data(), opts.stride,
						weights.posy.size() * sizeof(decltype(weights.posy)::value_type));
			flush_range(weights.posz.data(), opts.stride,
						weights.posz.size() * sizeof(decltype(weights.posz)::value_type));
	}

	void animate_joints()
	{
		animate_joints_scattered_ownerless_soa(root);
	}

	void animate_weights()
	{
		animate_weights_scattered_parent_soa(weights);
	}

	void reset_weights()
	{
		for (size_t i = 0; i < weights.w.size(); i++) {
			weights.posx[i] = initial[i].x;
			weights.posy[i] = initial[i].y;
			weights.posz[i] = initial[i].z;
		}
	}
};

struct PooledOwnerlessSoa
{
	std::vector<JointPooledOwnerlessSoa> joints;
	WeightPooledParentSoa weights;
	std::vector<vec3> initial;

	void initialize_from_opts(BenchmarkOptions& opts)
	{
		joints.resize(opts.joint_count);

		weights.w.reserve(opts.weight_count);
		weights.posx.reserve(opts.weight_count);
		weights.posy.reserve(opts.weight_count);
		weights.posz.reserve(opts.weight_count);
		weights.joint.reserve(opts.weight_count);

		auto& root = joints[0];
		root.orient_relative = quat::from_xyz(
			opts.joint_data[0].qx, opts.joint_data[0].qy, opts.joint_data[0].qz);

		root.orient_absolute.w = 1;
		root.orient_absolute.x = 0;
		root.orient_absolute.y = 0;
		root.orient_absolute.z = 0;

		for (size_t i = 1; i < opts.joint_count; i++) {
			auto& new_joint = joints[i];
			new_joint.orient_relative = quat::from_xyz(
				opts.joint_data[i].qx, opts.joint_data[i].qy, opts.joint_data[i].qz);

			new_joint.orient_absolute.w = 1;
			new_joint.orient_absolute.x = 0;
			new_joint.orient_absolute.y = 0;
			new_joint.orient_absolute.z = 0;

			new_joint.parent = &joints[opts.joint_data[i].parent];
		}

		for (size_t c = 0; c < opts.weight_replication_count; c++) {
			for (size_t i = 0; i < opts.weight_count; i++) {
				auto* joint = &joints[opts.weight_data[i].joint];

				weights.w.push_back(opts.weight_data[i].w);
				weights.posx.push_back(opts.weight_data[i].x);
				weights.posy.push_back(opts.weight_data[i].y);
				weights.posz.push_back(opts.weight_data[i].z);
				weights.joint.push_back(joint);

				vec3 v;
				v.x = opts.weight_data[i].x;
				v.y = opts.weight_data[i].y;
				v.z = opts.weight_data[i].z;

				initial.push_back(v);
			}
		}
	}

	void flush_from_cache(const BenchmarkOptions& opts)
	{
		flush_range(joints.data(), opts.stride,
					joints.size() * sizeof(decltype(joints)::value_type));

		flush_range(weights.joint.data(), opts.stride,
					weights.joint.size() * sizeof(decltype(weights.joint)::value_type));
		flush_range(weights.w.data(), opts.stride,
					weights.w.size() * sizeof(decltype(weights.w)::value_type));
		flush_range(weights.posx.data(), opts.stride,
					weights.posx.size() * sizeof(decltype(weights.posx)::value_type));
		flush_range(weights.posy.data(), opts.stride,
					weights.posy.size() * sizeof(decltype(weights.posy)::value_type));
		flush_range(weights.posz.data(), opts.stride,
					weights.posz.size() * sizeof(decltype(weights.posz)::value_type));
	}

	void animate_joints()
	{
		animate_joints_pooled_ownerless_soa(joints);
	}

	void animate_weights()
	{
		animate_weights_pooled_parent_soa(weights);
	}

	void reset_weights()
	{
		for (size_t i = 0; i < weights.w.size(); i++) {
			weights.posx[i] = initial[i].x;
			weights.posy[i] = initial[i].y;
			weights.posz[i] = initial[i].z;
		}
	}
};

template<typename DataStructure>
void BM_Template(benchmark::State& state)
{
	auto opts = from_args(state);
	DataStructure structure;
	structure.initialize_from_opts(opts);

	structure.animate_joints();
	structure.flush_from_cache(opts);

	for (auto _: state) {
		structure.animate_weights();
		structure.flush_from_cache(opts);
		benchmark::ClobberMemory();

		state.PauseTiming();
		structure.reset_weights();
		structure.flush_from_cache(opts);
		state.ResumeTiming();
	}

	state.SetComplexityN(opts.weight_count * opts.weight_replication_count);
	auto items = state.iterations() * opts.weight_count * opts.weight_replication_count;
	state.SetItemsProcessed(items);
	state.SetBytesProcessed(3 * sizeof(float) * items);
}

double max_of_vector(const std::vector<double>& v)
{
	return *(std::max_element(v.begin(), v.end()));
}

double min_of_vector(const std::vector<double>& v)
{
	return *(std::min_element(v.begin(), v.end()));
}

BENCHMARK_TEMPLATE(BM_Template, ScatteredOwningAos)
	->ArgNames({"Replication", "Dataset"})
	->Apply(CustomArgs)
	->Complexity(benchmark::oN)
	->ComputeStatistics("min", min_of_vector)
	->ComputeStatistics("max", max_of_vector);

BENCHMARK_TEMPLATE(BM_Template, PooledOwningAos)
	->ArgNames({"Replication", "Dataset"})
	->Apply(CustomArgs)
	->Complexity(benchmark::oN)
	->ComputeStatistics("min", min_of_vector)
	->ComputeStatistics("max", max_of_vector);

BENCHMARK_TEMPLATE(BM_Template, ScatteredOwnerlessAos)
	->ArgNames({"Replication", "Dataset"})
	->Apply(CustomArgs)
	->Complexity(benchmark::oN)
	->ComputeStatistics("min", min_of_vector)
	->ComputeStatistics("max", max_of_vector);

BENCHMARK_TEMPLATE(BM_Template, PooledOwnerlessAos)
	->ArgNames({"Replication", "Dataset"})
	->Apply(CustomArgs)
	->Complexity(benchmark::oN)
	->ComputeStatistics("min", min_of_vector)
	->ComputeStatistics("max", max_of_vector);

BENCHMARK_TEMPLATE(BM_Template, ScatteredOwningSoa)
	->ArgNames({"Replication", "Dataset"})
	->Apply(CustomArgs)
	->Complexity(benchmark::oN)
	->ComputeStatistics("min", min_of_vector)
	->ComputeStatistics("max", max_of_vector);

BENCHMARK_TEMPLATE(BM_Template, PooledOwningSoa)
	->ArgNames({"Replication", "Dataset"})
	->Apply(CustomArgs)
	->Complexity(benchmark::oN)
	->ComputeStatistics("min", min_of_vector)
	->ComputeStatistics("max", max_of_vector);

BENCHMARK_TEMPLATE(BM_Template, ScatteredOwnerlessSoa)
	->ArgNames({"Replication", "Dataset"})
	->Apply(CustomArgs)
	->Complexity(benchmark::oN)
	->ComputeStatistics("min", min_of_vector)
	->ComputeStatistics("max", max_of_vector);

BENCHMARK_TEMPLATE(BM_Template, PooledOwnerlessSoa)
	->ArgNames({"Replication", "Dataset"})
	->Apply(CustomArgs)
	->Complexity(benchmark::oN)
	->ComputeStatistics("min", min_of_vector)
	->ComputeStatistics("max", max_of_vector);

BENCHMARK_MAIN();
