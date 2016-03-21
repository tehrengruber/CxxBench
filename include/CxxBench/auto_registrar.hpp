namespace CxxBench {

struct AutoRegistrar {
	AutoRegistrar(std::string name, bench_case_t bench_case) {
		CxxBench::registry.add(name, bench_case);
	}
};

} // end CxxBench