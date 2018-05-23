#include "anim_smart.h"

#include <benchmark/benchmark.h>
#include <vector>

static void animate_smart_joint_impl(benchmark::State& state, JointSmart& joint)
{
	for (auto& e: joint.children) {
		e->orient_absolute = e->orient_relative * joint.orient_absolute;
	}

	for (auto& e: joint.children) {
		animate_smart_joint_impl(state, *e);
	}
}

void animate_smart(benchmark::State& state, JointSmart& joint, std::vector<WeightSmart>& weights)
{
	state.PauseTiming();

	joint.orient_absolute = joint.orient_relative;
	animate_smart_joint_impl(state, joint);

	state.ResumeTiming();

	for (auto& e: weights) {
		auto* joint = e.joint;
		e.pos = joint->orient_absolute * e.pos;
	}
}
