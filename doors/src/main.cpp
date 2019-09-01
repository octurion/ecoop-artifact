#include <array>
#include <cstdint>
#include <random>
#include <vector>

#include <benchmark/benchmark.h>

static constexpr size_t NUM_DOORS = 100;
static constexpr size_t NUM_CHARACTERS = 100;

static constexpr float DIMENSION_MAX = 100;
static constexpr float DOOR_DIST = 1.0f;

static const std::array<uint32_t, 100> SEEDS {
	14159, 26535, 89793, 23846, 26433, 83279, 50288, 41971, 69399, 37510,
	58209, 74944, 59230, 78164,  6286, 20899, 86280, 34825, 34211, 70679,
	82148,  8651, 32823,  6647,  9384, 46095, 50582, 23172, 53594,  8128,
	48111, 74502, 84102, 70193, 85211,  5559, 64462, 29489, 54930, 38196,
	44288, 10975, 66593, 34461, 28475, 64823, 37867, 83165, 27120, 19091,
	45648, 56692, 34603, 48610, 45432, 66482, 13393, 60726,  2491, 41273,
	72458, 70066,  6315, 58817, 48815, 20920, 96282, 92540, 91715, 36436,
	78925, 90360,  1133,  5305, 48820, 46652, 13841, 46951, 94151, 16094,
	33057, 27036, 57595, 91953,  9218, 61173, 81932, 61179, 31051, 18548,
	 7446, 23799, 62749, 56735, 18857, 52724, 89122, 79381, 83011, 94912,
};

enum class Allegiance: uint8_t
{
	Red, Blue
};

struct DoorOnePool
{
	float x = 0;
	float y = 0;
	Allegiance allegiance = Allegiance::Red;
	bool open = false;
};

struct DoorManyPools
{
	float x = 0;
	float y = 0;

	Allegiance allegiance = Allegiance::Red;
	bool open = false;
};

struct CharacterOnePool
{
	float x = 0;
	float y = 0;

	Allegiance allegiance = Allegiance::Red;
};

struct CharacterManyPools
{
	float x = 0;
	float y = 0;

	Allegiance allegiance = Allegiance::Red;
};

void reset_doors(std::vector<DoorOnePool>& pool)
{
	for (auto& e: pool) {
		e.open = false;
	}
}

void reset_doors(std::vector<DoorManyPools>& pool)
{
	for (auto& e: pool) {
		e.open = false;
	}
}

void open_doors(std::vector<DoorOnePool>& doors, std::vector<CharacterOnePool>& characters)
{
	for (auto& d: doors) {
		for (auto& c: characters) {
			if (c.allegiance != d.allegiance) {
				continue;
			}

			float dx = c.x - d.x;
			float dy = c.y - d.y;

			float dist2 = dx * dx + dy * dy;
			if (dist2 > DOOR_DIST * DOOR_DIST) {
				continue;
			}

			d.open = true;
			break;
		}
	}
}

void open_doors(std::vector<DoorManyPools>& doors, std::vector<CharacterManyPools>& characters)
{
	for (auto& d: doors) {
		for (auto& c: characters) {
			float dx = c.x - d.x;
			float dy = c.y - d.y;

			float dist2 = dx * dx + dy * dy;
			if (dist2 > DOOR_DIST * DOOR_DIST) {
				continue;
			}

			d.open = true;
			break;
		}
	}
}

size_t num_open_doors(const std::vector<DoorOnePool>& doors)
{
	size_t count = 0;
	for (const auto& e: doors) {
		if (e.open) {
			count++;
		}
	}

	return count;
}

size_t num_open_doors(const std::vector<DoorManyPools>& doors)
{
	size_t count = 0;
	for (const auto& e: doors) {
		if (e.open) {
			count++;
		}
	}

	return count;
}

void generate_entities(
	std::vector<DoorOnePool>& doors,
	std::vector<CharacterOnePool>& characters,
	uint32_t seed,
	float allegiance_probability)
{
	std::mt19937 mt(seed);
	std::uniform_real_distribution<float> coord_dist(0, DIMENSION_MAX);
	std::uniform_real_distribution<float> allegiance_dist(0, 1);

	for (size_t i = 0; i < NUM_DOORS; i++) {
		DoorOnePool door;
		door.x = coord_dist(mt);
		door.y = coord_dist(mt);
		door.allegiance = (Allegiance) (allegiance_dist(mt) < allegiance_probability);
		door.open = false;

		doors.push_back(door);
	}

	for (size_t i = 0; i < NUM_CHARACTERS; i++) {
		CharacterOnePool character;
		character.x = coord_dist(mt);
		character.y = coord_dist(mt);
		character.allegiance = (Allegiance) (allegiance_dist(mt) < allegiance_probability);

		characters.push_back(character);
	}
}

void generate_entities(
	std::vector<DoorManyPools>& red_doors,
	std::vector<DoorManyPools>& blue_doors,
	std::vector<CharacterManyPools>& red_characters,
	std::vector<CharacterManyPools>& blue_characters,
	uint64_t seed,
	float allegiance_probability)
{
	std::mt19937 mt(seed);
	std::uniform_real_distribution<float> coord_dist(0, DIMENSION_MAX);
	std::uniform_real_distribution<float> allegiance_dist(0, 1);

	for (size_t i = 0; i < NUM_DOORS; i++) {
		DoorManyPools door;
		door.x = coord_dist(mt);
		door.y = coord_dist(mt);
		door.open = false;
		door.allegiance = (Allegiance) (allegiance_dist(mt) < allegiance_probability);

		if (door.allegiance == Allegiance::Red) {
			red_doors.push_back(door);
		} else {
			blue_doors.push_back(door);
		}
	}

	for (size_t i = 0; i < NUM_CHARACTERS; i++) {
		CharacterManyPools character;
		character.x = coord_dist(mt);
		character.y = coord_dist(mt);
		character.allegiance = (Allegiance) (allegiance_dist(mt) < allegiance_probability);

		if (character.allegiance == Allegiance::Red) {
			red_characters.push_back(character);
		} else {
			blue_characters.push_back(character);
		}
	}
}

void BM_DoorsOnePool(benchmark::State& state)
{
	uint64_t seed = state.range(0);
	float allegiance_probability = state.range(1) / 100.f;

	std::vector<DoorOnePool> doors;
	std::vector<CharacterOnePool> characters;

	generate_entities(doors, characters, seed, allegiance_probability);

	for (auto _: state) {
		open_doors(doors, characters);
		benchmark::DoNotOptimize(num_open_doors(doors));

		state.PauseTiming();
		reset_doors(doors);
		state.ResumeTiming();
	}
}

void BM_DoorsManyPools(benchmark::State& state)
{
	uint64_t seed = state.range(0);
	float allegiance_probability = state.range(1) / 100.f;

	std::vector<DoorManyPools> red_doors;
	std::vector<DoorManyPools> blue_doors;
	std::vector<CharacterManyPools> red_characters;
	std::vector<CharacterManyPools> blue_characters;

	generate_entities(red_doors, blue_doors, red_characters, blue_characters, seed, allegiance_probability);

	for (auto _: state) {
		open_doors(red_doors, red_characters);
		open_doors(blue_doors, blue_characters);

		benchmark::DoNotOptimize(num_open_doors(red_doors) + num_open_doors(blue_doors));

		state.PauseTiming();
		reset_doors(red_doors);
		reset_doors(blue_doors);
		state.ResumeTiming();
	}
}

static void CustomArguments(benchmark::internal::Benchmark* b) {
	for (const auto& j: SEEDS) {
		b->Args({(uint32_t) j, 50});
		b->Args({(uint32_t) j, 70});
		b->Args({(uint32_t) j, 90});
	}
}

BENCHMARK(BM_DoorsOnePool)
	->Apply(CustomArguments)
	->ArgNames({"Mersenne seed", "Allegiance probability"});

BENCHMARK(BM_DoorsManyPools)
	->Apply(CustomArguments)
	->ArgNames({"Mersenne seed", "Allegiance probability"});

BENCHMARK_MAIN();
