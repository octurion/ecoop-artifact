#pragma once

#include "anim.h"

#include <benchmark/benchmark.h>
#include <memory>
#include <vector>

struct JointSmart
{
	quat orient_relative;
	quat orient_absolute;

	std::vector<std::unique_ptr<JointSmart>> children;
};

struct WeightSmart
{
	JointSmart* joint;
	vec3 pos;
	float w;
};

void animate_smart(benchmark::State& state, JointSmart& joint, std::vector<WeightSmart>& weights);
