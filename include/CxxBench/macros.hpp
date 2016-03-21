#define CXX_MACRO_HELPER_CAT(prefix, suffix)	prefix ## suffix
#define _CXX_BENCH_UNIQUE_NAME(prefix, suffix)	CXX_MACRO_HELPER_CAT(prefix, suffix)
#define CXX_BENCH_UNIQUE_NAME(prefix)			_CXX_BENCH_UNIQUE_NAME(cxx_bench_case_ ## prefix, __LINE__)

#define CXX_BENCHMARK_CASE( BenchCaseName ) \
	static void BenchCaseName(CxxBench::State& state); \
	CxxBench::AutoRegistrar CXX_BENCH_UNIQUE_NAME(BenchCaseName) (#BenchCaseName, &BenchCaseName); \
	static void BenchCaseName(CxxBench::State& state)