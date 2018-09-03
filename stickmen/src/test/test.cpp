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
	bool flush_cache;
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

	opt.flush_cache = state.range(2);

	return opt;
}

static void CustomArgs(benchmark::internal::Benchmark* b)
{
	for (int i = 1; i <= 10; i++) {
			b->Args({i * 32, 0, true});
		for (int j = 0; j < static_cast<int>(DatasetType::DatasetCount); j++) {
			// b->Args({i, j, true});
		}
	}
}

static void flush_cache(const void* ptr, size_t stride, size_t len)
{
	size_t flush_count = (len + stride) & (stride - 1);

	for (size_t i = 0; i < flush_count; i += stride) {
		_mm_clflush((const char*) ptr + i * stride);
	}
}

static void flush_cache_scattered_owning(JointScatteredOwning& root, size_t stride)
{
	for (auto& e: root.children) {
		flush_cache_scattered_owning(*e, stride);
	}

	flush_cache(root.weights.data(),
				stride,
				root.weights.size() * sizeof(decltype(root.weights)::value_type));

	flush_cache(root.children.data(),
				stride,
				root.children.size() * sizeof(decltype(root.children)::value_type));

	flush_cache(&root, stride, sizeof(root));
}

static void BM_ScatteredOwning(benchmark::State& state) {
	size_t stride = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);
	std::unordered_map<int, JointScatteredOwning*> joint_map;

	JointScatteredOwning root;
	joint_map[0] = &root;

	auto opts = from_args(state);

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

	for (auto _: state) {
		state.PauseTiming();
		animate_joints_scattered_owning(root);
		state.ResumeTiming();

		animate_weights_scattered_parentless(root);
		benchmark::ClobberMemory();

		state.PauseTiming();
		if (opts.flush_cache) {
			flush_cache_scattered_owning(root, stride);
		}
		state.ResumeTiming();

	}

	state.SetItemsProcessed(state.iterations() * opts.weight_count * opts.weight_replication_count);
}
BENCHMARK(BM_ScatteredOwning)
	->ArgNames({"Replication", "Dataset", "FlushCache"})
	->Apply(CustomArgs);

static void flush_cache_scattered_owning_soa(JointScatteredOwningSoa& root, size_t stride)
{
	for (auto& e: root.children) {
		flush_cache_scattered_owning_soa(*e, stride);
	}

	auto& weights = root.weights;

	flush_cache(weights.w.data(),
				stride,
				weights.w.size() * sizeof(decltype(weights.w)::value_type));
	flush_cache(weights.posx.data(),
				stride,
				weights.posx.size() * sizeof(decltype(weights.posx)::value_type));
	flush_cache(weights.posy.data(),
				stride,
				weights.posy.size() * sizeof(decltype(weights.posy)::value_type));
	flush_cache(weights.posz.data(),
				stride,
				weights.posz.size() * sizeof(decltype(weights.posz)::value_type));

	flush_cache(root.children.data(),
				stride,
				root.children.size() * sizeof(decltype(root.children)::value_type));

	flush_cache(&root, stride, sizeof(root));
}

static void BM_ScatteredOwningSoa(benchmark::State& state) {
	size_t stride = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);
	std::unordered_map<int, JointScatteredOwningSoa*> joint_map;

	JointScatteredOwningSoa root;
	joint_map[0] = &root;

	auto opts = from_args(state);

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
		}
	}

	for (auto _: state) {
		state.PauseTiming();
		animate_joints_scattered_owning_soa(root);
		state.ResumeTiming();

		animate_weights_scattered_parentless_soa(root);
		benchmark::ClobberMemory();

		state.PauseTiming();
		if (opts.flush_cache) {
			flush_cache_scattered_owning_soa(root, stride);
		}
		state.ResumeTiming();

	}
	state.SetItemsProcessed(state.iterations() * opts.weight_count * opts.weight_replication_count);
}
BENCHMARK(BM_ScatteredOwningSoa)
	->ArgNames({"Replication", "Dataset", "FlushCache"})
	->Apply(CustomArgs);

static void flush_cache_scattered_ownerless_joints(JointScatteredOwnerless& root, size_t stride)
{
	for (auto& e: root.children) {
		flush_cache_scattered_ownerless_joints(*e, stride);
	}

	flush_cache(root.children.data(),
				stride,
				root.children.size() * sizeof(decltype(root.children)::value_type));

	flush_cache(&root, stride, sizeof(root));
}

static void BM_ScatteredOwnerless(benchmark::State& state) {
	size_t stride = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);
	std::unordered_map<int, JointScatteredOwnerless*> joint_map;

	JointScatteredOwnerless root;
	joint_map[0] = &root;

	auto opts = from_args(state);
	std::vector<WeightScatteredParent> weights;
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
		}
	}

	for (auto _: state) {
		state.PauseTiming();
		animate_joints_scattered_ownerless(root);
		state.ResumeTiming();

		animate_weights_scattered_parent(weights);
		benchmark::ClobberMemory();

		state.PauseTiming();
		if (opts.flush_cache) {
			flush_cache_scattered_ownerless_joints(root, stride);
			flush_cache(weights.data(), stride,
						weights.size() * sizeof(decltype(weights)::value_type));
		}
		state.ResumeTiming();
	}
	state.SetItemsProcessed(state.iterations() * opts.weight_count * opts.weight_replication_count);
}
BENCHMARK(BM_ScatteredOwnerless)
	->ArgNames({"Replication", "Dataset", "FlushCache"})
	->Apply(CustomArgs);

static void flush_cache_scattered_ownerless_soa_joints(JointScatteredOwnerlessSoa& root, size_t stride)
{
	for (auto& e: root.children) {
		flush_cache_scattered_ownerless_soa_joints(*e, stride);
	}

	flush_cache(root.children.data(),
				stride,
				root.children.size() * sizeof(decltype(root.children)::value_type));

	flush_cache(&root, stride, sizeof(root));
}
static void BM_ScatteredOwnerlessSoa(benchmark::State& state) {
	size_t stride = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);
	std::unordered_map<int, JointScatteredOwnerlessSoa*> joint_map;

	JointScatteredOwnerlessSoa root;
	joint_map[0] = &root;

	auto opts = from_args(state);

	WeightScatteredParentSoa weights;
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
		}
	}

	for (auto _: state) {
		state.PauseTiming();
		animate_joints_scattered_ownerless_soa(root);
		state.ResumeTiming();

		animate_weights_scattered_parent_soa(weights);
		benchmark::ClobberMemory();

		state.PauseTiming();
		if (opts.flush_cache) {
			flush_cache_scattered_ownerless_soa_joints(root, stride);
			flush_cache(weights.joint.data(), stride,
						weights.joint.size() * sizeof(decltype(weights.joint)::value_type));
			flush_cache(weights.w.data(), stride,
						weights.w.size() * sizeof(decltype(weights.w)::value_type));
			flush_cache(weights.posx.data(), stride,
						weights.posx.size() * sizeof(decltype(weights.posx)::value_type));
			flush_cache(weights.posy.data(), stride,
						weights.posy.size() * sizeof(decltype(weights.posy)::value_type));
			flush_cache(weights.posz.data(), stride,
						weights.posz.size() * sizeof(decltype(weights.posz)::value_type));
		}
		state.ResumeTiming();
	}
	state.SetItemsProcessed(state.iterations() * opts.weight_count * opts.weight_replication_count);
}
BENCHMARK(BM_ScatteredOwnerlessSoa)
	->ArgNames({"Replication", "Dataset", "FlushCache"})
	->Apply(CustomArgs);

static void BM_PooledOwning(benchmark::State& state) {
	size_t stride = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);

	auto opts = from_args(state);

	std::vector<JointPooledOwning> joints(opts.joint_count);

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

	for (auto _: state) {
		state.PauseTiming();
		animate_joints_pooled_owning(joints);
		state.ResumeTiming();

		animate_weights_pooled_parentless(joints);
		benchmark::ClobberMemory();

		state.PauseTiming();
		if (opts.flush_cache) {
			for (auto& e: joints) {
				flush_cache(e.weights.data(), stride,
							e.weights.size() * sizeof(decltype(e.weights)::value_type));
			}
			flush_cache(joints.data(), stride,
						joints.size() * sizeof(decltype(joints)::value_type));
		}
		state.ResumeTiming();
	}
	state.SetItemsProcessed(state.iterations() * opts.weight_count * opts.weight_replication_count);
}
BENCHMARK(BM_PooledOwning)
	->ArgNames({"Replication", "Dataset", "FlushCache"})
	->Apply(CustomArgs);

static void BM_PooledOwningSoa(benchmark::State& state) {
	size_t stride = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);

	auto opts = from_args(state);

	std::vector<JointPooledOwningSoa> joints(opts.joint_count);

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
		}
	}

	for (auto _: state) {
		state.PauseTiming();
		animate_joints_pooled_owning_soa(joints);
		state.ResumeTiming();

		animate_weights_pooled_parentless_soa(joints);
		benchmark::ClobberMemory();

		state.PauseTiming();
		if (opts.flush_cache) {
			for (auto& e: joints) {
				flush_cache(e.weights.w.data(), stride,
							e.weights.w.size() * sizeof(decltype(e.weights.w)::value_type));

				flush_cache(e.weights.posx.data(), stride,
							e.weights.posx.size() * sizeof(decltype(e.weights.posx)::value_type));
				flush_cache(e.weights.posy.data(), stride,
							e.weights.posy.size() * sizeof(decltype(e.weights.posy)::value_type));
				flush_cache(e.weights.posz.data(), stride,
							e.weights.posz.size() * sizeof(decltype(e.weights.posz)::value_type));
			}
			flush_cache(joints.data(), stride,
						joints.size() * sizeof(decltype(joints)::value_type));
		}
		state.ResumeTiming();
	}
	state.SetItemsProcessed(state.iterations() * opts.weight_count * opts.weight_replication_count);
}
BENCHMARK(BM_PooledOwningSoa)
	->ArgNames({"Replication", "Dataset", "FlushCache"})
	->Apply(CustomArgs);

static void BM_PooledOwnerless(benchmark::State& state) {
	size_t stride = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);

	auto opts = from_args(state);

	std::vector<JointPooledOwnerless> joints(opts.joint_count);
	std::vector<WeightPooledParent> weights;
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
		}
	}

	for (auto _: state) {
		state.PauseTiming();
		animate_joints_pooled_ownerless(joints);
		state.ResumeTiming();

		animate_weights_pooled_parent(weights);
		benchmark::ClobberMemory();

		state.PauseTiming();
		if (opts.flush_cache) {
			flush_cache(joints.data(), stride,
						joints.size() * sizeof(decltype(joints)::value_type));
			flush_cache(weights.data(), stride,
						weights.size() * sizeof(decltype(weights)::value_type));
		}
		state.ResumeTiming();
	}
	state.SetItemsProcessed(state.iterations() * opts.weight_count * opts.weight_replication_count);
}
BENCHMARK(BM_PooledOwnerless)
	->ArgNames({"Replication", "Dataset", "FlushCache"})
	->Apply(CustomArgs);

static void BM_PooledOwnerlessSoa(benchmark::State& state) {
	size_t stride = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);

	auto opts = from_args(state);

	std::vector<JointPooledOwnerlessSoa> joints(opts.joint_count);

	WeightPooledParentSoa weights;
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
		}
	}

	for (auto _: state) {
		state.PauseTiming();
		animate_joints_pooled_ownerless_soa(joints);
		state.ResumeTiming();

		animate_weights_pooled_parent_soa(weights);
		benchmark::ClobberMemory();

		state.PauseTiming();
		if (opts.flush_cache) {
			flush_cache(joints.data(), stride,
						joints.size() * sizeof(decltype(joints)::value_type));

			flush_cache(weights.joint.data(), stride,
						weights.joint.size() * sizeof(decltype(weights.joint)::value_type));
			flush_cache(weights.w.data(), stride,
						weights.w.size() * sizeof(decltype(weights.w)::value_type));
			flush_cache(weights.posx.data(), stride,
						weights.posx.size() * sizeof(decltype(weights.posx)::value_type));
			flush_cache(weights.posy.data(), stride,
						weights.posy.size() * sizeof(decltype(weights.posy)::value_type));
			flush_cache(weights.posz.data(), stride,
						weights.posz.size() * sizeof(decltype(weights.posz)::value_type));
		}
		state.ResumeTiming();
	}
	state.SetItemsProcessed(state.iterations() * opts.weight_count * opts.weight_replication_count);
}
BENCHMARK(BM_PooledOwnerlessSoa)
	->ArgNames({"Replication", "Dataset", "FlushCache"})
	->Apply(CustomArgs);

static void BM_PooledOwnerlessScattered(benchmark::State& state) {
	size_t stride = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);

	auto opts = from_args(state);

	std::vector<JointPooledOwnerlessScattered> joints(opts.joint_count);
	std::vector<std::unique_ptr<WeightPooledParentScattered>> weights;

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
			std::unique_ptr<WeightPooledParentScattered> new_weight(new WeightPooledParentScattered);

			new_weight->w = opts.weight_data[i].w;
			new_weight->pos.x = opts.weight_data[i].x;
			new_weight->pos.y = opts.weight_data[i].y;
			new_weight->pos.z = opts.weight_data[i].z;
			new_weight->joint = joint;

			weights.emplace_back(std::move(new_weight));
		}
	}

	for (auto _: state) {
		state.PauseTiming();
		animate_joints_pooled_ownerless_scattered(joints);
		state.ResumeTiming();

		animate_weights_pooled_parent_scattered(weights);
		benchmark::ClobberMemory();

		state.PauseTiming();
		if (opts.flush_cache) {
			flush_cache(joints.data(), stride,
						joints.size() * sizeof(decltype(joints)::value_type));

			for (auto& e: weights) {
				flush_cache (e.get(), stride, sizeof(*e));
			}

			flush_cache(weights.data(), stride,
						weights.size() * sizeof(decltype(weights)::value_type));
		}
		state.ResumeTiming();
	}
	state.SetItemsProcessed(state.iterations() * opts.weight_count * opts.weight_replication_count);
}
BENCHMARK(BM_PooledOwnerlessScattered)
	->ArgNames({"Replication", "Dataset", "FlushCache"})
	->Apply(CustomArgs);

static void BM_PooledOwnerScattered(benchmark::State& state) {
	size_t stride = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);

	auto opts = from_args(state);

	std::vector<JointPooledOwningScattered> joints(opts.joint_count);

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
			std::unique_ptr<WeightParentless> new_weight(new WeightParentless);

			new_weight->w = opts.weight_data[i].w;
			new_weight->pos.x = opts.weight_data[i].x;
			new_weight->pos.y = opts.weight_data[i].y;
			new_weight->pos.z = opts.weight_data[i].z;

			joint->weights.emplace_back(std::move(new_weight));
		}
	}

	for (auto _: state) {
		state.PauseTiming();
		animate_joints_pooled_owning_scattered(joints);
		state.ResumeTiming();

		animate_weights_pooled_parentless_scattered(joints);
		benchmark::ClobberMemory();

		state.PauseTiming();
		if (opts.flush_cache) {
			for (auto& e: joints) {
				for (auto& w: e.weights) {
					flush_cache (w.get(), stride, sizeof(*w));
				}
				flush_cache(e.weights.data(), stride,
							e.weights.size() * sizeof(decltype(e.weights)::value_type));
			}
			flush_cache(joints.data(), stride,
						joints.size() * sizeof(decltype(joints)::value_type));
		}
		state.ResumeTiming();
	}
	state.SetItemsProcessed(state.iterations() * opts.weight_count * opts.weight_replication_count);
}
BENCHMARK(BM_PooledOwnerScattered)
	->ArgNames({"Replication", "Dataset", "FlushCache"})
	->Apply(CustomArgs);

BENCHMARK_MAIN();
