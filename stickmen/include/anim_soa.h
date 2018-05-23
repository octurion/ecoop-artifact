#pragma once

#include "anim.h"

#include <benchmark/benchmark.h>
#include <vector>

struct JointSoa
{
	quat orient_relative;
	quat orient_absolute;

	std::vector<JointSoa> children;
};

struct WeightSoa
{
	std::vector<JointSoa*> joint;

	std::vector<float> posx;
	std::vector<float> posy;
	std::vector<float> posz;

	std::vector<float> w;
};

void animate_soa(benchmark::State& state, JointSoa& joint, WeightSoa& weights);
