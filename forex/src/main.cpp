#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

#include <benchmark/benchmark.h>

#define xstr(s) str(s)
#define str(s) #s

#define RATE_MAX_LEN 15
#define ISO_8601_LEN 10
#define INPUT_FILE "../eurofxref-hist.csv"

#define DATE_RECENT "2018-01-01"
#define PROB_RECENT 0.8

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

struct Rate
{
	char date[ISO_8601_LEN + 1];
	double USD, JPY, BGN, CYP, CZK, DKK, EEK, GBP, HUF, LTL, LVL, MTL, PLN,
		   ROL, RON, SEK, SIT, SKK, CHF, ISK, NOK, HRK, RUB, TRL, TRY, AUD,
		   BRL, CAD, CNY, HKD, IDR, ILS, INR, KRW, MXN, MYR, NZD, PHP, SGD,
		   THB, ZAR;
};

struct RateFrequent
{
	char date[ISO_8601_LEN + 1];
	double USD, GBP;
};
struct RateInfrequent
{
	double JPY, BGN, CYP, CZK, DKK, EEK, HUF, LTL, LVL, MTL, PLN, ROL, RON,
		   SEK, SIT, SKK, CHF, ISK, NOK, HRK, RUB, TRL, TRY, AUD, BRL, CAD,
		   CNY, HKD, IDR, ILS, INR, KRW, MXN, MYR, NZD, PHP, SGD, THB, ZAR;
};

enum class Currency : uint8_t
{
	USD, GBP
};

struct Query
{
	char date[ISO_8601_LEN + 1];
	Currency currency;
};

template<typename Rate>
struct Comp
{
	bool operator()(const Rate& rate, const char* date) const {
		return strcmp(date, rate.date) < 0;
	}
};

static bool initialized = false;
static std::vector<Rate> rates;

static void read_rate(FILE* in, double* ptr)
{
	char rate_str[RATE_MAX_LEN + 1];

	int unused = fscanf(in, "%" xstr(RATE_MAX_LEN) "[^,],", rate_str);
	(void) unused;

	errno = 0;
	char* marker = rate_str;
	double val = strtod(rate_str, &marker);
	if (errno != 0 || marker == rate_str) {
		val = HUGE_VAL;
	}

	*ptr = val;
}

static void read_exchange_rates()
{
	if (initialized) {
		return;
	}

	// Ensure that '.' is the decimal separator used
	setlocale(LC_ALL, "C");

	FILE* in = fopen(INPUT_FILE, "r");
	if (in == NULL) {
		perror(INPUT_FILE);
		abort();
	}

	// First line is the header, skip it
	int unused = fscanf(in, "%*s\n");
	(void) unused;

	for (;;) {
		Rate rate;

		if (fscanf(in, "%" xstr(ISO_8601_LEN) "[^,],", rate.date) == EOF) {
			break;
		}

		read_rate(in, &rate.USD);
		read_rate(in, &rate.JPY);
		read_rate(in, &rate.BGN);
		read_rate(in, &rate.CYP);
		read_rate(in, &rate.CZK);
		read_rate(in, &rate.DKK);
		read_rate(in, &rate.EEK);
		read_rate(in, &rate.GBP);
		read_rate(in, &rate.HUF);
		read_rate(in, &rate.LTL);
		read_rate(in, &rate.LVL);
		read_rate(in, &rate.MTL);
		read_rate(in, &rate.PLN);
		read_rate(in, &rate.ROL);
		read_rate(in, &rate.RON);
		read_rate(in, &rate.SEK);
		read_rate(in, &rate.SIT);
		read_rate(in, &rate.SKK);
		read_rate(in, &rate.CHF);
		read_rate(in, &rate.ISK);
		read_rate(in, &rate.NOK);
		read_rate(in, &rate.HRK);
		read_rate(in, &rate.RUB);
		read_rate(in, &rate.TRL);
		read_rate(in, &rate.TRY);
		read_rate(in, &rate.AUD);
		read_rate(in, &rate.BRL);
		read_rate(in, &rate.CAD);
		read_rate(in, &rate.CNY);
		read_rate(in, &rate.HKD);
		read_rate(in, &rate.IDR);
		read_rate(in, &rate.ILS);
		read_rate(in, &rate.INR);
		read_rate(in, &rate.KRW);
		read_rate(in, &rate.MXN);
		read_rate(in, &rate.MYR);
		read_rate(in, &rate.NZD);
		read_rate(in, &rate.PHP);
		read_rate(in, &rate.SGD);
		read_rate(in, &rate.THB);
		read_rate(in, &rate.ZAR);

		rates.push_back(rate);

		unused = fscanf(in, "\n");
		(void) unused;
	}

	std::reverse(rates.begin(), rates.end());

	fclose(in);

	initialized = true;
}

class Fixture : public benchmark::Fixture
{
public:
	std::vector<Query> queries;

	void SetUp(const ::benchmark::State& state) override {
		read_exchange_rates();

		auto from_2018_onward = std::lower_bound(
			rates.begin(), rates.end(), DATE_RECENT,
			[](const Rate& rate, const char* date) {
				return strcmp(date, rate.date) < 0;
			});

		size_t size_old = std::distance(rates.begin(), from_2018_onward);
		size_t size_new = std::distance(from_2018_onward, rates.end());

		std::mt19937 rng(state.range(1));
		std::uniform_int_distribution<int> fifty_fifty(0, 1);
		std::uniform_real_distribution<double> prob_dist(0, 1);

		std::uniform_int_distribution<size_t> old_dist(0, size_old + 1);
		std::uniform_int_distribution<size_t> new_dist(0, size_new + 1);

		size_t range = state.range();
		for (size_t i = 0; i < range; i++) {
			Query q;

			bool recent_date = prob_dist(rng) < PROB_RECENT;
			auto it = recent_date
				? rates.begin() + old_dist(rng)
				: from_2018_onward + new_dist(rng);

			strcpy(q.date, it->date);

			bool usd = fifty_fifty(rng);
			q.currency = usd ? Currency::USD : Currency::GBP;

			queries.push_back(q);
		}
	}

	void TearDown(const ::benchmark::State&) override {
		queries.clear();
	}
};

BENCHMARK_DEFINE_F(Fixture, AosOnePool)(benchmark::State& st)
{
	for (auto _ : st) {
		for (const auto& q: queries) {
			auto it = std::lower_bound(rates.begin(), rates.end(), q.date, Comp<Rate>());

			if (it == rates.end()) {
				continue;
			}

			auto rate = q.currency == Currency::USD
				? it->USD
				: it->GBP;

			benchmark::DoNotOptimize(rate);
		}
	}
	st.SetComplexityN(st.range());
}

BENCHMARK_DEFINE_F(Fixture, MixedManyPools)(benchmark::State& st)
{
	auto from_2018_onward = std::lower_bound(
		rates.begin(), rates.end(), DATE_RECENT, Comp<Rate>());

	std::vector<RateFrequent> frequent;
	std::vector<RateInfrequent> infrequent;

	std::vector<Rate> rare;

	for (auto it = from_2018_onward; it != rates.end(); it++) {
		RateFrequent freq;
		strcpy(freq.date, it->date);
		freq.USD = it->USD;
		freq.GBP = it->GBP;

		frequent.push_back(freq);

		RateInfrequent infreq;
		infreq.JPY = it->JPY;
		infreq.BGN = it->BGN;
		infreq.CYP = it->CYP;
		infreq.CZK = it->CZK;
		infreq.DKK = it->DKK;
		infreq.EEK = it->EEK;
		infreq.HUF = it->HUF;
		infreq.LTL = it->LTL;
		infreq.LVL = it->LVL;
		infreq.MTL = it->MTL;
		infreq.PLN = it->PLN;
		infreq.ROL = it->ROL;
		infreq.RON = it->RON;
		infreq.SEK = it->SEK;
		infreq.SIT = it->SIT;
		infreq.SKK = it->SKK;
		infreq.CHF = it->CHF;
		infreq.ISK = it->ISK;
		infreq.NOK = it->NOK;
		infreq.HRK = it->HRK;
		infreq.RUB = it->RUB;
		infreq.TRL = it->TRL;
		infreq.TRY = it->TRY;
		infreq.AUD = it->AUD;
		infreq.BRL = it->BRL;
		infreq.CAD = it->CAD;
		infreq.CNY = it->CNY;
		infreq.HKD = it->HKD;
		infreq.IDR = it->IDR;
		infreq.ILS = it->ILS;
		infreq.INR = it->INR;
		infreq.KRW = it->KRW;
		infreq.MXN = it->MXN;
		infreq.MYR = it->MYR;
		infreq.NZD = it->NZD;
		infreq.PHP = it->PHP;
		infreq.SGD = it->SGD;
		infreq.THB = it->THB;
		infreq.ZAR = it->ZAR;

		infrequent.push_back(infreq);
	}
	for (auto it = rates.begin(); it != from_2018_onward; it++) {
		rare.push_back(*it);
	}

	for (auto _ : st) {
		for (const auto& q: queries) {
			if (strcmp(q.date, DATE_RECENT) >= 0) {
				auto freq_it = std::lower_bound(
					frequent.begin(), frequent.end(), q.date, Comp<RateFrequent>());
				if (freq_it != frequent.end()) {
					auto rate = q.currency == Currency::USD
						? freq_it->USD
						: freq_it->GBP;

					benchmark::DoNotOptimize(rate);
				}
			} else {
				auto rare_it = std::lower_bound(
					rare.begin(), rare.end(), q.date, Comp<Rate>());
				if (rare_it != rare.end()) {
					auto rate = q.currency == Currency::USD
						? rare_it->USD
						: rare_it->GBP;

					benchmark::DoNotOptimize(rate);
				}
			}
		}
	}
	st.SetComplexityN(st.range());
}

static void CustomArguments(benchmark::internal::Benchmark* b) {
	for (const auto& j: SEEDS) {
		for (size_t i = 0; i < 10; i++) {
			b->Args({8 << i, (int64_t) j});
		}
	}
}

BENCHMARK_REGISTER_F(Fixture, AosOnePool)
	->Apply(CustomArguments)
	->Complexity(benchmark::oN);
BENCHMARK_REGISTER_F(Fixture, MixedManyPools)
	->Apply(CustomArguments)
	->Complexity(benchmark::oN);

BENCHMARK_MAIN();
