#include <array>
#include <random>
#include <vector>

#include <benchmark/benchmark.h>

static constexpr size_t NUM_DOORS = 100;
static constexpr size_t NUM_CHARACTERS = 100;

static constexpr float DIMENSION_MAX = 100;
static constexpr float DOOR_DIST = 1.0f;

static const std::array<uint64_t, 100> SEEDS {
	1415926535ULL, 8979323846ULL, 2643383279ULL, 5028841971ULL, 6939937510ULL,
	5820974944ULL, 5923078164ULL,  628620899ULL, 8628034825ULL, 3421170679ULL,
	8214808651ULL, 3282306647ULL,  938446095ULL, 5058223172ULL, 5359408128ULL,
	4811174502ULL, 8410270193ULL, 8521105559ULL, 6446229489ULL, 5493038196ULL,
	4428810975ULL, 6659334461ULL, 2847564823ULL, 3786783165ULL, 2712019091ULL,
	4564856692ULL, 3460348610ULL, 4543266482ULL, 1339360726ULL,  249141273ULL,
	7245870066ULL,  631558817ULL, 4881520920ULL, 9628292540ULL, 9171536436ULL,
	7892590360ULL,  113305305ULL, 4882046652ULL, 1384146951ULL, 9415116094ULL,
	3305727036ULL, 5759591953ULL,  921861173ULL, 8193261179ULL, 3105118548ULL,
	 744623799ULL, 6274956735ULL, 1885752724ULL, 8912279381ULL, 8301194912ULL,
	9833673362ULL, 4406566430ULL, 8602139494ULL, 6395224737ULL, 1907021798ULL,
	6094370277ULL,  539217176ULL, 2931767523ULL, 8467481846ULL, 7669405132ULL,
	   5681271ULL, 4526356082ULL, 7785771342ULL, 7577896091ULL, 7363717872ULL,
	1468440901ULL, 2249534301ULL, 4654958537ULL, 1050792279ULL, 6892589235ULL,
	4201995611ULL, 2129021960ULL, 8640344181ULL, 5981362977ULL, 4771309960ULL,
	5187072113ULL, 4999999837ULL, 2978049951ULL,  597317328ULL, 1609631859ULL,
	5024459455ULL, 3469083026ULL, 4252230825ULL, 3344685035ULL, 2619311881ULL,
	7101000313ULL, 7838752886ULL, 5875332083ULL, 8142061717ULL, 7669147303ULL,
	5982534904ULL, 2875546873ULL, 1159562863ULL, 8823537875ULL, 9375195778ULL,
	1857780532ULL, 1712268066ULL, 1300192787ULL, 6611195909ULL, 2164201989ULL,
};

enum class Allegiance
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
			if (dist2 <= DOOR_DIST * DOOR_DIST) {
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
			if (dist2 <= DOOR_DIST * DOOR_DIST) {
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

void generate_entities(std::vector<DoorOnePool>& doors, std::vector<CharacterOnePool>& characters, uint64_t seed)
{
	std::mt19937 mt(seed);
	std::uniform_real_distribution<float> coord_dist(0, DIMENSION_MAX);
	std::uniform_int_distribution<int> allegiance_dist(0, 1);

	for (size_t i = 0; i < NUM_DOORS; i++) {
		DoorOnePool door;
		door.x = coord_dist(mt);
		door.y = coord_dist(mt);
		door.allegiance = (Allegiance) allegiance_dist(mt);
		door.open = false;

		doors.push_back(door);
	}

	for (size_t i = 0; i < NUM_CHARACTERS; i++) {
		CharacterOnePool character;
		character.x = coord_dist(mt);
		character.y = coord_dist(mt);
		character.allegiance = (Allegiance) allegiance_dist(mt);

		characters.push_back(character);
	}
}

void generate_entities(
	std::vector<DoorManyPools>& red_doors,
	std::vector<DoorManyPools>& blue_doors,
	std::vector<CharacterManyPools>& red_characters,
	std::vector<CharacterManyPools>& blue_characters,
	uint64_t seed)
{
	std::mt19937 mt(seed);
	std::uniform_real_distribution<float> coord_dist(0, DIMENSION_MAX);
	std::uniform_int_distribution<int> allegiance_dist(0, 1);

	for (size_t i = 0; i < NUM_DOORS; i++) {
		DoorManyPools door;
		door.x = coord_dist(mt);
		door.y = coord_dist(mt);
		door.open = false;

		if (allegiance_dist(mt) == 0) {
			red_doors.push_back(door);
		} else {
			blue_doors.push_back(door);
		}
	}

	for (size_t i = 0; i < NUM_CHARACTERS; i++) {
		CharacterManyPools character;
		character.x = coord_dist(mt);
		character.y = coord_dist(mt);

		if (allegiance_dist(mt) == 0) {
			red_characters.push_back(character);
		} else {
			blue_characters.push_back(character);
		}
	}
}

void BM_DoorsOnePool(benchmark::State& state)
{
	uint64_t seed = state.range(0);

	std::vector<DoorOnePool> doors;
	std::vector<CharacterOnePool> characters;

	generate_entities(doors, characters, seed);

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

	std::vector<DoorManyPools> red_doors;
	std::vector<DoorManyPools> blue_doors;
	std::vector<CharacterManyPools> red_characters;
	std::vector<CharacterManyPools> blue_characters;

	generate_entities(red_doors, blue_doors, red_characters, blue_characters, seed);

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
		b->Args({(int64_t) j});
	}
}

BENCHMARK(BM_DoorsOnePool)
	->Apply(CustomArguments)
	->ArgNames({"Mersenne seed"});

BENCHMARK(BM_DoorsManyPools)
	->Apply(CustomArguments)
	->ArgNames({"Mersenne seed"});

BENCHMARK_MAIN();
