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
#define INPUT_FILE "forex/eurofxref-hist.csv"

#define DATE_RECENT "2018-01-01"
#define PROB_RECENT 0.8

#define NUM_QUERIES 5000

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

struct DateComp
{
	bool operator()(const std::array<char, ISO_8601_LEN + 1>& rate, const char* date) const {
		return strcmp(rate.data(), date) < 0;
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

BENCHMARK_DEFINE_F(Fixture, MixedManyPoolsSoa)(benchmark::State& st)
{
	auto from_2018_onward = std::lower_bound(
		rates.begin(), rates.end(), DATE_RECENT, Comp<Rate>());

	std::vector<std::array<char, ISO_8601_LEN + 1>> dates;
	std::vector<double> USD;
	std::vector<double> JPY;
	std::vector<double> BGN;
	std::vector<double> CYP;
	std::vector<double> CZK;
	std::vector<double> DKK;
	std::vector<double> EEK;
	std::vector<double> GBP;
	std::vector<double> HUF;
	std::vector<double> LTL;
	std::vector<double> LVL;
	std::vector<double> MTL;
	std::vector<double> PLN;
	std::vector<double> ROL;
	std::vector<double> RON;
	std::vector<double> SEK;
	std::vector<double> SIT;
	std::vector<double> SKK;
	std::vector<double> CHF;
	std::vector<double> ISK;
	std::vector<double> NOK;
	std::vector<double> HRK;
	std::vector<double> RUB;
	std::vector<double> TRL;
	std::vector<double> TRY;
	std::vector<double> AUD;
	std::vector<double> BRL;
	std::vector<double> CAD;
	std::vector<double> CNY;
	std::vector<double> HKD;
	std::vector<double> IDR;
	std::vector<double> ILS;
	std::vector<double> INR;
	std::vector<double> KRW;
	std::vector<double> MXN;
	std::vector<double> MYR;
	std::vector<double> NZD;
	std::vector<double> PHP;
	std::vector<double> SGD;
	std::vector<double> THB;
	std::vector<double> ZAR;

	std::vector<Rate> rare;

	for (auto it = from_2018_onward; it != rates.end(); it++) {
		std::array<char, ISO_8601_LEN + 1> date = {};
		strcpy(date.data(), it->date);
		dates.push_back(date);

		USD.push_back(it->USD);
		GBP.push_back(it->GBP);

		JPY.push_back(it->JPY);
		BGN.push_back(it->BGN);
		CYP.push_back(it->CYP);
		CZK.push_back(it->CZK);
		DKK.push_back(it->DKK);
		EEK.push_back(it->EEK);
		HUF.push_back(it->HUF);
		LTL.push_back(it->LTL);
		LVL.push_back(it->LVL);
		MTL.push_back(it->MTL);
		PLN.push_back(it->PLN);
		ROL.push_back(it->ROL);
		RON.push_back(it->RON);
		SEK.push_back(it->SEK);
		SIT.push_back(it->SIT);
		SKK.push_back(it->SKK);
		CHF.push_back(it->CHF);
		ISK.push_back(it->ISK);
		NOK.push_back(it->NOK);
		HRK.push_back(it->HRK);
		RUB.push_back(it->RUB);
		TRL.push_back(it->TRL);
		TRY.push_back(it->TRY);
		AUD.push_back(it->AUD);
		BRL.push_back(it->BRL);
		CAD.push_back(it->CAD);
		CNY.push_back(it->CNY);
		HKD.push_back(it->HKD);
		IDR.push_back(it->IDR);
		ILS.push_back(it->ILS);
		INR.push_back(it->INR);
		KRW.push_back(it->KRW);
		MXN.push_back(it->MXN);
		MYR.push_back(it->MYR);
		NZD.push_back(it->NZD);
		PHP.push_back(it->PHP);
		SGD.push_back(it->SGD);
		THB.push_back(it->THB);
		ZAR.push_back(it->ZAR);
	}
	for (auto it = rates.begin(); it != from_2018_onward; it++) {
		rare.push_back(*it);
	}

	for (auto _ : st) {
		for (const auto& q: queries) {
			if (strcmp(q.date, DATE_RECENT) >= 0) {
				auto freq_it = std::lower_bound(
					dates.begin(), dates.end(), q.date, DateComp());
				if (freq_it != dates.end()) {
					auto idx = std::distance(dates.begin(), freq_it);
					auto rate = q.currency == Currency::USD
						? USD[idx]
						: GBP[idx];

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
		b->Args({NUM_QUERIES, (uint32_t) j});
	}
}

BENCHMARK_REGISTER_F(Fixture, AosOnePool)
	->Apply(CustomArguments)
	->Complexity(benchmark::oN);
BENCHMARK_REGISTER_F(Fixture, MixedManyPools)
	->Apply(CustomArguments)
	->Complexity(benchmark::oN);
BENCHMARK_REGISTER_F(Fixture, MixedManyPoolsSoa)
	->Apply(CustomArguments)
	->Complexity(benchmark::oN);

BENCHMARK_MAIN();
