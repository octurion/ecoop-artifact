#pragma once

#include "anim.h"

#include <benchmark/benchmark.h>
#include <memory>
#include <vector>

struct WeightNaive
{
	vec3 pos;
	float w;
};

struct JointNaive
{
	quat orient_relative;
	quat orient_absolute;

	std::vector<std::unique_ptr<JointNaive>> children;
	std::vector<WeightNaive> weights;
};

void animate_naive(benchmark::State& state, JointNaive& joint);
