#pragma once

#define NUM_JOINTS1 110
#define NUM_WEIGHTS1 3206

#define NUM_JOINTS2 33
#define NUM_WEIGHTS2 867

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

extern const struct JointData JOINTS1[NUM_JOINTS1];
extern const struct WeightData WEIGHTS1[NUM_WEIGHTS1];

extern const struct JointData JOINTS2[NUM_JOINTS2];
extern const struct WeightData WEIGHTS2[NUM_WEIGHTS2];
