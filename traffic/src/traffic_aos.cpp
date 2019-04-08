#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdint>
#include <random>
#include <vector>

#include <benchmark/benchmark.h>

static const int MAX_VELOCITY = 10;
static const int MAX_DEGREE = 4;

static const int NUM_INTERSECTIONS = 20;
static const float CELL_LENGTH = 0.005f;
static const float PRODUCER_RATIO = 0.02f;
static const float TARGET_RATIO = 0.002f;
static const float CAR_ALLOCATION_RATIO = 0.02f;
static const float CELL_TARGET_RATIO = 0.002f;
static const float SLOW_DOWN_PROBABILITY = 0.2f;

static const int TRAFFIC_LIGHT_PHASE = 5;

static const size_t NULLPTR_IDX = -1;

static const std::array<uint64_t, 20> SEEDS {
	31415926535ULL, 8979323846ULL, 2643383279ULL, 5028841971ULL, 6939937510ULL,
	 5820974944ULL, 5923078164ULL,  628620899ULL, 8628034825ULL, 3421170679ULL,
	 8214808651ULL, 3282306647ULL,  938446095ULL, 5058223172ULL, 5359408128ULL,
	 4811174502ULL, 8410270193ULL, 8521105559ULL, 6446229489ULL, 5493038196ULL,
};

std::uniform_int_distribution<uint32_t> max_vel_dist(MAX_VELOCITY / 2, MAX_VELOCITY);
std::uniform_real_distribution<float> probability_dist(0.0, 1.0);
std::uniform_real_distribution<float> pos_dist(0.0, 1.0);

enum class CellFlags: uint32_t {
	Producer     = 1 << 0,
	HasCar       = 1 << 1,
	IsTarget     = 1 << 2,
	ShouldOccupy = 1 << 3,
};

struct Cell
{
	std::mt19937 rand_state;
	std::mt19937 car_rand_state;

	std::array<size_t, MAX_DEGREE> incoming;
	size_t num_incoming;

	std::array<size_t, MAX_DEGREE> outgoing;
	size_t num_outgoing;

	std::array<size_t, MAX_VELOCITY> car_path;
	size_t num_car_path;

	uint32_t car_velocity;
	uint32_t car_max_velocity;

	uint32_t curr_max_velocity;
	uint32_t max_velocity;

	float x;
	float y;

	uint32_t flags;

	template<typename Rng>
	Cell(uint32_t max_velocity, float x, float y, uint32_t flags, Rng& rng)
		: rand_state(rng())
		, car_rand_state(rng())
		, incoming()
		, num_incoming(0)
		, outgoing()
		, num_outgoing(0)
		, car_path()
		, num_car_path(0)
		, car_velocity(0)
		, car_max_velocity(0)
		, curr_max_velocity(0)
		, max_velocity(max_velocity)
		, x(x)
		, y(y)
		, flags(flags)
	{
	}
};

struct TrafficLight
{
	std::array<size_t, 2 * MAX_DEGREE> cells = {};
	size_t num_cells = 0;

	uint32_t timer = 0;
	uint32_t phase_time = 5;
	uint32_t phase = 0;
};

void create_street_network(std::vector<Cell>& cells, std::vector<TrafficLight>& traffic_lights, uint64_t seed)
{
	std::random_device rd;
	std::mt19937 mt(seed);

	for (size_t i = 0; i < NUM_INTERSECTIONS; i++) {
		float x = pos_dist(mt);
		float y = pos_dist(mt);

		uint32_t max_velocity = max_vel_dist(mt);

		cells.emplace_back(max_velocity, x, y, 0, mt);
	}

	std::uniform_int_distribution<size_t> intersection_dist(0, NUM_INTERSECTIONS - 1);
	std::uniform_int_distribution<size_t> outgoing_dist(1, MAX_DEGREE);

	for (size_t i = 0; i < NUM_INTERSECTIONS; i++) {
		size_t degree = outgoing_dist(mt);
		uint32_t max_velocity = cells[i].max_velocity;
		for (size_t idx = 0; idx < degree; idx++) {
			size_t j;
			do {
				j = intersection_dist(mt);
			} while (j == i || cells[j].num_incoming >= MAX_DEGREE);

			float dx = cells[i].x - cells[j].x;
			float dy = cells[i].y - cells[j].y;
			float dist = sqrt(dx * dx + dy * dy);
			size_t steps = dist / CELL_LENGTH;
			float step_x = dx / steps;
			float step_y = dy / steps;

			size_t last = i;
			for (size_t k = 0; k < steps; k++) {
				float new_x = cells[i].x + k * step_x;
				float new_y = cells[i].y + k * step_y;

				bool is_producer = probability_dist(mt) < PRODUCER_RATIO;
				uint32_t flags = is_producer ? (uint32_t) CellFlags::Producer : 0;

				bool is_target = probability_dist(mt) < CELL_TARGET_RATIO;
				flags &= (is_target ? (uint32_t) CellFlags::IsTarget : 0);

				cells.emplace_back(max_velocity, new_x, new_y, flags, mt);

				auto& old_cell = cells[last];
				old_cell.outgoing[old_cell.num_outgoing++] = cells.size() - 1;

				auto& new_cell = cells.back();
				new_cell.incoming[new_cell.num_incoming++] = last;

				last = cells.size() - 1;
			}

			auto& old_cell = cells[last];
			old_cell.outgoing[old_cell.num_outgoing++] = j;

			auto& new_cell = cells[j];
			new_cell.incoming[new_cell.num_incoming++] = last;
		}
	}

	for (size_t i = 0; i < NUM_INTERSECTIONS; i++) {
		traffic_lights.emplace_back();
		for (size_t j = 0; j < cells[i].num_incoming; j++) {
			auto cell_id = cells[i].incoming[j];

			traffic_lights[i].cells[traffic_lights[i].num_cells++] = cell_id;
			cells[cell_id].curr_max_velocity = 0;
		}
	}
}

void create_cars(std::vector<Cell>& cells)
{
	for (auto& e: cells) {
		uint32_t MASK = ((uint32_t) CellFlags::Producer) & ~(uint32_t) CellFlags::HasCar;

		bool can_create_car = (e.flags & MASK) != 0;
		if (!can_create_car) {
			continue;
		}

		bool must_create = probability_dist(e.rand_state) < PRODUCER_RATIO;
		if (!must_create) {
			continue;
		}

		e.flags &= (uint32_t) CellFlags::HasCar;
		e.num_car_path = 0;
		e.car_velocity = 0;
		e.max_velocity = max_vel_dist(e.rand_state);
	}
}

void step_traffic_lights(std::vector<Cell>& cells, std::vector<TrafficLight>& traffic_lights)
{
	for (auto& e: traffic_lights) {
		if (e.num_cells == 0) {
			continue;
		}

		e.timer++;
		e.timer = (e.timer == e.phase_time) ? 0 : e.timer;

		if (e.timer == 0) {
			cells[e.cells[e.phase]].curr_max_velocity = 0;

			e.phase++;
			e.phase = (e.phase == e.num_cells) ? 0 : e.phase;

			cells[e.cells[e.phase]].curr_max_velocity = cells[e.cells[e.phase]].max_velocity;
		}
	}
}

void prepare_paths(std::vector<Cell>& cells)
{
	std::uniform_int_distribution<uint32_t> speedup_path_dist(1, 2);
	for (size_t idx = 0; idx < cells.size(); idx++) {
		auto& e = cells[idx];
		if ((e.flags & (uint32_t) CellFlags::HasCar) != 0) {
			continue;
		}

		e.num_car_path = 0;

		auto speedup = speedup_path_dist(e.car_rand_state);
		e.car_velocity = std::max(e.car_max_velocity, e.car_velocity + speedup);

		{
			size_t curr_cell = idx;
			size_t next_cell = idx;
			size_t length = e.car_velocity;
			for (size_t i = 0; i < length; i++) {
				if ((cells[curr_cell].flags & (uint32_t) CellFlags::IsTarget) || cells[curr_cell].num_outgoing == 0) {
					break;
				}

				if (cells[curr_cell].num_outgoing > 1) {
					std::uniform_int_distribution<size_t> path_dist(0, cells[curr_cell].num_outgoing - 1);
					auto choice = path_dist(e.car_rand_state);
					next_cell = cells[curr_cell].outgoing[choice];
				} else {
					next_cell = cells[curr_cell].outgoing[0];
				}

				if ((cells[next_cell].flags & (uint32_t) CellFlags::HasCar) != 0) {
					break;
				}

				curr_cell = next_cell;
				e.car_path[e.num_car_path++] = curr_cell;
			}

			e.car_velocity = e.num_car_path;
		}

		{
			e.car_velocity = std::max(e.car_velocity, e.curr_max_velocity);
			size_t path_index = 0;
			size_t distance = 1;

			while (distance <= e.car_velocity) {
				size_t next_cell = e.car_path[path_index];

				if ((cells[next_cell].flags & (uint32_t) CellFlags::HasCar) != 0) {
					distance--;
					e.car_velocity = distance;
					break;
				}

				if (e.car_velocity > cells[next_cell].curr_max_velocity) {
					if (cells[next_cell].curr_max_velocity > distance - 1) {
						e.car_velocity = cells[next_cell].curr_max_velocity;
					} else {
						distance--;
						e.car_velocity = distance;
						break;
					}
				}

				distance++;
				path_index++;
			}
			distance--;
		}

		if (probability_dist(e.car_rand_state) < SLOW_DOWN_PROBABILITY) {
			e.car_velocity--;
		}
	}
}

void step_move(std::vector<Cell>& cells)
{
	for (auto& e: cells) {
		if ((e.flags & (uint32_t) CellFlags::HasCar) == 0) {
			continue;
		}
		if (e.car_velocity == 0) {
			continue;
		}

		auto& occupied = cells[e.car_path[e.num_car_path - 1]];
		occupied.flags &= (uint32_t) CellFlags::ShouldOccupy;
		occupied.car_velocity = e.car_velocity;
		occupied.car_max_velocity = e.car_max_velocity;
		occupied.car_rand_state = e.car_rand_state;

		std::copy(e.car_path.begin(), e.car_path.begin() + e.num_car_path,
				  occupied.car_path.begin());
		occupied.num_car_path = e.num_car_path;

		e.flags ^= (uint32_t) CellFlags::HasCar;
	}
}

void commit_occupy(std::vector<Cell>& cells)
{
	for (auto& e: cells) {
		if ((e.flags & (uint32_t) CellFlags::ShouldOccupy) != 0) {
			e.flags ^= (uint32_t) CellFlags::ShouldOccupy;
			e.flags &= (uint32_t) CellFlags::HasCar;
		}

		if (e.num_outgoing == 0 || (e.flags & (uint32_t) CellFlags::IsTarget) != 0) {
			e.flags &= ~((uint32_t) CellFlags::HasCar);
		}
	}
}

void BM_TrafficAos(benchmark::State& state)
{
	for (auto _: state) {
		state.PauseTiming();

		std::vector<Cell> cells;
		std::vector<TrafficLight> traffic_lights;

		create_street_network(cells, traffic_lights, (uint64_t) state.range(1));

		size_t num_iterations = state.range(0);
		state.ResumeTiming();

		for (size_t i = 0; i < num_iterations; i++) {
			create_cars(cells);
			step_traffic_lights(cells, traffic_lights);
			prepare_paths(cells);
			step_move(cells);
			commit_occupy(cells);
		}
	}
}

static void CustomArguments(benchmark::internal::Benchmark* b) {
	for (const auto& j: SEEDS) {
		b->Args({100, (int64_t)j});
		b->Args({500, (int64_t)j});
		b->Args({1000, (int64_t)j});
	}
}

BENCHMARK(BM_TrafficAos)
	->Apply(CustomArguments)
	->ArgNames({"Num iterations", "Mersenne seed"});

BENCHMARK_MAIN();
