#include "anim_naive.h"

#include <benchmark/benchmark.h>
#include <vector>

static void animate_naive_impl(benchmark::State& state, JointNaive& joint)
{
	for (auto& e: joint.children) {
		e.orient_absolute = e.orient_relative * joint.orient_absolute;
	}

	for (auto& e: joint.children) {
		animate_naive_impl(state, e);
	}

	state.ResumeTiming();

	for (auto& e: joint.weights) {
		e.pos = joint.orient_absolute * e.pos;
	}

	state.PauseTiming();
}

void animate_naive(benchmark::State& state, JointNaive& joint)
{
	state.PauseTiming();

	joint.orient_absolute = joint.orient_relative;
	animate_naive_impl(state, joint);

	state.ResumeTiming();
}
