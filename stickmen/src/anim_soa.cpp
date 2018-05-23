#include "anim_soa.h"

#include <benchmark/benchmark.h>
#include <vector>

static void animate_soa_joint_impl(benchmark::State& state, JointSoa& joint)
{
	for (auto& e: joint.children) {
		e->orient_absolute = e->orient_relative * joint.orient_absolute;
	}

	for (auto& e: joint.children) {
		animate_soa_joint_impl(state, *e);
	}
}

void animate_soa(benchmark::State& state, JointSoa& joint, WeightSoa& weights)
{
	state.PauseTiming();

	joint.orient_absolute = joint.orient_relative;
	animate_soa_joint_impl(state, joint);

	state.ResumeTiming();

	size_t len = weights.joint.size();
	for (size_t i = 0; i < len; i++) {
		auto* joint = weights.joint[i];

		vec3 pos;
		pos.x = weights.posx[i];
		pos.y = weights.posy[i];
		pos.z = weights.posz[i];

		auto new_pos = joint->orient_absolute * pos;

		weights.posx[i] = new_pos.x;
		weights.posy[i] = new_pos.y;
		weights.posz[i] = new_pos.z;
	}
}
