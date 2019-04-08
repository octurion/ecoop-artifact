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
	31415926535UL, 8979323846UL, 2643383279UL, 5028841971UL, 6939937510UL,
	 5820974944UL, 5923078164UL,  628620899UL, 8628034825UL, 3421170679UL,
	 8214808651UL, 3282306647UL,  938446095UL, 5058223172UL, 5359408128UL,
	 4811174502UL, 8410270193UL, 8521105559UL, 6446229489UL, 5493038196UL,
};

std::uniform_int_distribution<uint32_t> max_vel_dist(MAX_VELOCITY / 2, MAX_VELOCITY);
std::uniform_real_distribution<float> probability_dist(0.0, 1.0);
std::uniform_real_distribution<float> pos_dist(0.0, 1.0);

struct Xoroshiro128State
{
	uint64_t state0;
	uint64_t state1;
};

template<typename Rng>
void generate_xoroshiro_state(Xoroshiro128State& state, Rng& rng)
{
	std::uniform_int_distribution<uint64_t> u64_distr(0, UINT64_MAX);
	state.state0 = u64_distr(rng);
	state.state1 = u64_distr(rng);
}

enum class CellFlags: uint32_t {
	Producer     = 1 << 0,
	HasCar       = 1 << 1,
	IsTarget     = 1 << 2,
	ShouldOccupy = 1 << 3,
};

template<typename T, size_t N>
struct VecArray
{
	std::array<T, N> elems = {};
	size_t size = 0;
};

struct CellPool
{
	std::vector<std::mt19937> rand_state;
	std::vector<std::mt19937> car_rand_state;

	std::vector<VecArray<size_t, MAX_DEGREE>> incoming;
	std::vector<VecArray<size_t, MAX_DEGREE>> outgoing;
	std::vector<VecArray<size_t, MAX_VELOCITY>> car_path;

	std::vector<uint32_t> car_velocity;
	std::vector<uint32_t> car_max_velocity;

	std::vector<uint32_t> curr_max_velocity;
	std::vector<uint32_t> max_velocity;

	std::vector<float> x;
	std::vector<float> y;

	std::vector<uint32_t> flags;

	template<typename Rng>
	void addCell(uint32_t max_velocity, float x, float y, uint32_t flags, Rng& rng)
	{
		this->rand_state.emplace_back(rng());
		this->car_rand_state.emplace_back(rng());
		this->incoming.emplace_back();
		this->outgoing.emplace_back();
		this->car_path.emplace_back();
		this->car_velocity.emplace_back(0);
		this->car_max_velocity.emplace_back(0);
		this->curr_max_velocity.emplace_back(0);
		this->max_velocity.emplace_back(max_velocity);
		this->x.emplace_back(x);
		this->y.emplace_back(y);
		this->flags.emplace_back(flags);
	}
};

struct TrafficLightPool
{
	std::vector<VecArray<size_t, 2 * MAX_DEGREE>> cells;

	std::vector<uint32_t> timer;
	std::vector<uint32_t> phase_time;
	std::vector<uint32_t> phase;

	void addTrafficLight() {
		cells.emplace_back();
		timer.emplace_back(0);
		phase_time.emplace_back(5);
		phase.emplace_back(0);
	}
};

void create_street_network(CellPool& cells, TrafficLightPool& traffic_lights, uint64_t seed)
{
	std::random_device rd;
	std::mt19937 mt(seed);

	for (size_t i = 0; i < NUM_INTERSECTIONS; i++) {
		float x = pos_dist(mt);
		float y = pos_dist(mt);

		uint32_t max_velocity = max_vel_dist(mt);

		cells.addCell(max_velocity, x, y, 0, mt);
	}

	std::uniform_int_distribution<size_t> intersection_dist(0, NUM_INTERSECTIONS - 1);
	std::uniform_int_distribution<size_t> outgoing_dist(1, MAX_DEGREE);

	for (size_t i = 0; i < NUM_INTERSECTIONS; i++) {
		size_t degree = outgoing_dist(mt);
		uint32_t max_velocity = cells.max_velocity[i];
		for (size_t idx = 0; idx < degree; idx++) {
			size_t j;
			do {
				j = intersection_dist(mt);
			} while (j == i || cells.incoming[j].size >= MAX_DEGREE);

			float dx = cells.x[i] - cells.x[j];
			float dy = cells.y[i] - cells.y[j];
			float dist = sqrt(dx * dx + dy * dy);
			size_t steps = dist / CELL_LENGTH;
			float step_x = dx / steps;
			float step_y = dy / steps;

			size_t last = i;
			for (size_t k = 0; k < steps; k++) {
				float new_x = cells.x[i] + k * step_x;
				float new_y = cells.y[i] + k * step_y;

				bool is_producer = probability_dist(mt) < PRODUCER_RATIO;
				uint32_t flags = is_producer ? (uint32_t) CellFlags::Producer : 0;

				bool is_target = probability_dist(mt) < CELL_TARGET_RATIO;
				flags &= (is_target ? (uint32_t) CellFlags::IsTarget : 0);

				cells.addCell(max_velocity, new_x, new_y, flags, mt);

				auto& old_cell = cells.outgoing[last];
				old_cell.elems[old_cell.size++] = cells.x.size() - 1;

				auto& new_cell = cells.incoming.back();
				new_cell.elems[new_cell.size++] = last;

				last = cells.x.size() - 1;
			}

			auto& old_cell = cells.outgoing[last];
			old_cell.elems[old_cell.size++] = j;

			auto& new_cell = cells.incoming[j];
			new_cell.elems[new_cell.size++] = last;
		}
	}

	for (size_t i = 0; i < NUM_INTERSECTIONS; i++) {
		traffic_lights.addTrafficLight();
		for (size_t j = 0; j < cells.incoming[i].size; j++) {
			auto cell_id = cells.incoming[i].elems[j];

			traffic_lights.cells[i].elems[traffic_lights.cells[i].size++] = cell_id;
			cells.curr_max_velocity[cell_id] = 0;
		}
	}
}

void create_cars(CellPool& cells)
{
	for (size_t i = 0; i < cells.x.size(); i++) {
		uint32_t MASK = ((uint32_t) CellFlags::Producer) & ~(uint32_t) CellFlags::HasCar;

		bool can_create_car = (cells.flags[i] & MASK) != 0;
		if (!can_create_car) {
			continue;
		}

		bool must_create = probability_dist(cells.rand_state[i]) < PRODUCER_RATIO;
		if (!must_create) {
			continue;
		}

		cells.flags[i] &= (uint32_t) CellFlags::HasCar;
		cells.car_path[i].size = 0;
		cells.car_velocity[i] = 0;
		cells.max_velocity[i] = max_vel_dist(cells.rand_state[i]);
	}
}

void step_traffic_lights(CellPool& cells, TrafficLightPool& traffic_lights)
{
	for (size_t i = 0; i < traffic_lights.phase.size(); i++) {
		if (traffic_lights.cells[i].size == 0) {
			continue;
		}

		traffic_lights.timer[i]++;
		traffic_lights.timer[i] = (traffic_lights.timer[i] == traffic_lights.phase_time[i])
			? 0
			: traffic_lights.timer[i];

		if (traffic_lights.timer[i] == 0) {
			cells.curr_max_velocity[traffic_lights.cells[i].elems[traffic_lights.phase[i]]] = 0;

			traffic_lights.phase[i]++;
			traffic_lights.phase[i] = (traffic_lights.phase[i] == traffic_lights.cells[i].size)
				? 0
				: traffic_lights.phase[i];

			cells.curr_max_velocity[traffic_lights.cells[i].elems[traffic_lights.phase[i]]] =
				cells.max_velocity[traffic_lights.cells[i].elems[traffic_lights.phase[i]]];
		}
	}
}

void prepare_paths(CellPool& cells)
{
	std::uniform_int_distribution<uint32_t> speedup_path_dist(1, 2);
	for (size_t idx = 0; idx < cells.x.size(); idx++) {
		if ((cells.flags[idx] & (uint32_t) CellFlags::HasCar) != 0) {
			continue;
		}

		cells.car_path[idx].size = 0;

		auto speedup = speedup_path_dist(cells.car_rand_state[idx]);
		cells.car_velocity[idx] = std::max(cells.car_max_velocity[idx], cells.car_velocity[idx] + speedup);

		{
			size_t curr_cell = idx;
			size_t next_cell = idx;
			size_t length = cells.car_velocity[idx];
			for (size_t i = 0; i < length; i++) {
				if ((cells.flags[curr_cell] & (uint32_t) CellFlags::IsTarget) || cells.outgoing[curr_cell].size == 0) {
					break;
				}

				if (cells.outgoing[curr_cell].size > 1) {
					std::uniform_int_distribution<size_t> path_dist(0, cells.outgoing[curr_cell].size - 1);
					auto choice = path_dist(cells.car_rand_state[idx]);
					next_cell = cells.outgoing[curr_cell].elems[choice];
				} else {
					next_cell = cells.outgoing[curr_cell].elems[0];
				}

				if ((cells.flags[next_cell] & (uint32_t) CellFlags::HasCar) != 0) {
					break;
				}

				curr_cell = next_cell;
				cells.car_path[idx].elems[cells.car_path[idx].size++] = curr_cell;
			}

			cells.car_velocity[idx] = cells.car_path[idx].size;
		}

		{
			cells.car_velocity[idx] = std::max(cells.car_velocity[idx], cells.curr_max_velocity[idx]);
			size_t path_index = 0;
			size_t distance = 1;

			while (distance <= cells.car_velocity[idx]) {
				size_t next_cell = cells.car_path[idx].elems[path_index];
				
				if ((cells.flags[next_cell] & (uint32_t) CellFlags::HasCar) != 0) {
					distance--;
					cells.car_velocity[idx] = distance;
					break;
				}

				if (cells.car_velocity[idx] > cells.curr_max_velocity[next_cell]) {
					if (cells.curr_max_velocity[next_cell] > distance - 1) {
						cells.car_velocity[idx] = cells.curr_max_velocity[next_cell];
					} else {
						distance--;
						cells.car_velocity[idx] = distance;
						break;
					}
				}

				distance++;
				path_index++;
			}
			distance--;
		}

		if (probability_dist(cells.car_rand_state[idx]) < SLOW_DOWN_PROBABILITY) {
			cells.car_velocity[idx]--;
		}
	}
}

void step_move(CellPool& cells)
{
	for (size_t i = 0; i < cells.x.size(); i++) {
		if ((cells.flags[i] & (uint32_t) CellFlags::HasCar) == 0) {
			continue;
		}
		if (cells.car_velocity[i] == 0) {
			continue;
		}

		size_t occupied = cells.car_path[i].elems[cells.car_path[i].size - 1];
		cells.flags[occupied] &= (uint32_t) CellFlags::ShouldOccupy;
		cells.car_velocity[occupied] = cells.car_velocity[i];
		cells.car_max_velocity[occupied] = cells.car_max_velocity[i];
		cells.car_rand_state[occupied] = cells.car_rand_state[i];

		std::copy(cells.car_path[i].elems.begin(), cells.car_path[i].elems.begin() + cells.car_path[i].size,
				  cells.car_path[occupied].elems.begin());
		cells.car_path[occupied].size = cells.car_path[i].size;

		cells.flags[i] ^= (uint32_t) CellFlags::HasCar;
	}
}

void commit_occupy(CellPool& cells)
{
	for (size_t i = 0; i < cells.x.size(); i++) {
		if ((cells.flags[i] & (uint32_t) CellFlags::ShouldOccupy) != 0) {
			cells.flags[i] ^= (uint32_t) CellFlags::ShouldOccupy;
			cells.flags[i] &= (uint32_t) CellFlags::HasCar;
		}

		if (cells.outgoing[i].size == 0 || (cells.flags[i] & (uint32_t) CellFlags::IsTarget) != 0) {
			cells.flags[i] &= ~((uint32_t) CellFlags::HasCar);
		}
	}
}

void BM_TrafficSoa(benchmark::State& state)
{
	for (auto _: state) {
		state.PauseTiming();

		CellPool cells;
		TrafficLightPool traffic_lights;

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
	for (size_t i = 5; i <= 500; i *= 10) {
		for (const auto& j: SEEDS) {
			b->Args({(int64_t)i, (int64_t)j});
		}
	}
}

BENCHMARK(BM_TrafficSoa)
	->Apply(CustomArguments)
	->ArgNames({"Num iterations", "Mersenne seed"});

BENCHMARK_MAIN();
