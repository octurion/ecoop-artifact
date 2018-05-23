#include "anim_naive.h"
#include "anim_smart.h"
#include "anim_soa.h"

#include "test.h"

#include <benchmark/benchmark.h>

#include <unordered_map>
#include <utility>

static void BM_AnimNaive(benchmark::State& state) {
	std::unordered_map<int, JointNaive*> joint_map;

	JointNaive root;
	joint_map[0] = &root;

	root.orient_relative = quat::from_xyz(
		JOINTS[0].qx, JOINTS[0].qy, JOINTS[0].qz);

	root.orient_absolute.w = 1;
	root.orient_absolute.x = 0;
	root.orient_absolute.y = 0;
	root.orient_absolute.z = 0;

	for (size_t i = 1; i < NUM_JOINTS; i++) {
		auto& parent = joint_map[JOINTS[i].parent];

		std::unique_ptr<JointNaive> new_joint(new JointNaive);
		new_joint->orient_relative = quat::from_xyz(
			JOINTS[i].qx, JOINTS[i].qy, JOINTS[i].qz);

		new_joint->orient_absolute.w = 1;
		new_joint->orient_absolute.x = 0;
		new_joint->orient_absolute.y = 0;
		new_joint->orient_absolute.z = 0;

		joint_map[i] = new_joint.get();
		parent->children.emplace_back(std::move(new_joint));
	}

	for (size_t i = 0; i < NUM_WEIGHTS; i++) {
		auto& joint = joint_map[WEIGHTS[i].joint];

		WeightNaive new_weight;

		new_weight.w = WEIGHTS[i].w;

		new_weight.pos.x = WEIGHTS[i].x;
		new_weight.pos.y = WEIGHTS[i].y;
		new_weight.pos.z = WEIGHTS[i].z;

		joint->weights.emplace_back(std::move(new_weight));
	}

	while (state.KeepRunning()) {
		animate_naive(state, root);
	}

}
BENCHMARK(BM_AnimNaive);

static void BM_AnimSmart(benchmark::State& state) {
	std::unordered_map<int, JointSmart*> joint_map;

	std::vector<WeightSmart> weights;

	JointSmart root;
	joint_map[0] = &root;

	root.orient_relative = quat::from_xyz(
		JOINTS[0].qx, JOINTS[0].qy, JOINTS[0].qz);

	root.orient_absolute.w = 1;
	root.orient_absolute.x = 0;
	root.orient_absolute.y = 0;
	root.orient_absolute.z = 0;

	for (size_t i = 1; i < NUM_JOINTS; i++) {
		auto& parent = joint_map[JOINTS[i].parent];

		std::unique_ptr<JointSmart> new_joint(new JointSmart);
		new_joint->orient_relative = quat::from_xyz(
			JOINTS[i].qx, JOINTS[i].qy, JOINTS[i].qz);

		new_joint->orient_absolute.w = 1;
		new_joint->orient_absolute.x = 0;
		new_joint->orient_absolute.y = 0;
		new_joint->orient_absolute.z = 0;

		joint_map[i] = new_joint.get();
		parent->children.emplace_back(std::move(new_joint));
	}

	for (size_t i = 0; i < NUM_WEIGHTS; i++) {
		auto& joint = joint_map[WEIGHTS[i].joint];

		WeightSmart new_weight;

		new_weight.joint = joint;
		new_weight.w = WEIGHTS[i].w;

		new_weight.pos.x = WEIGHTS[i].x;
		new_weight.pos.y = WEIGHTS[i].y;
		new_weight.pos.z = WEIGHTS[i].z;

		weights.emplace_back(std::move(new_weight));
	}

	while (state.KeepRunning()) {
		animate_smart(state, root, weights);
	}
}
BENCHMARK(BM_AnimSmart);

static void BM_AnimSoa(benchmark::State& state) {
	std::unordered_map<int, JointSoa*> joint_map;

	JointSoa root;
	joint_map[0] = &root;

	root.orient_relative = quat::from_xyz(
		JOINTS[0].qx, JOINTS[0].qy, JOINTS[0].qz);

	root.orient_absolute.w = 1;
	root.orient_absolute.x = 0;
	root.orient_absolute.y = 0;
	root.orient_absolute.z = 0;

	for (size_t i = 1; i < NUM_JOINTS; i++) {
		auto& parent = joint_map[JOINTS[i].parent];

		std::unique_ptr<JointSoa> new_joint(new JointSoa);
		new_joint->orient_relative = quat::from_xyz(
			JOINTS[i].qx, JOINTS[i].qy, JOINTS[i].qz);

		new_joint->orient_absolute.w = 1;
		new_joint->orient_absolute.x = 0;
		new_joint->orient_absolute.y = 0;
		new_joint->orient_absolute.z = 0;

		joint_map[i] = new_joint.get();
		parent->children.emplace_back(std::move(new_joint));
	}

	WeightSoa weights;
	for (size_t i = 0; i < NUM_WEIGHTS; i++) {
		auto& joint = joint_map[WEIGHTS[i].joint];

		weights.joint.emplace_back(joint);
		weights.w.emplace_back(WEIGHTS[i].w);

		weights.posx.emplace_back(WEIGHTS[i].x);
		weights.posy.emplace_back(WEIGHTS[i].y);
		weights.posz.emplace_back(WEIGHTS[i].z);
	}

	while (state.KeepRunning()) {
		animate_soa(state, root, weights);
	}
}
BENCHMARK(BM_AnimSoa);

BENCHMARK_MAIN();
