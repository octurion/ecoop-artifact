#pragma once

#define NUM_JOINTS 110
#define NUM_WEIGHTS 3206

struct JointData
{
	int parent;

	float qx;
	float qy;
	float qz;
};

struct WeightData
{
	int joint;

	float w;

	float x;
	float y;
	float z;
};

extern const struct JointData JOINTS[NUM_JOINTS];
extern const struct WeightData WEIGHTS[NUM_WEIGHTS];
