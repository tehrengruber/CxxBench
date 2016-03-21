#include <vector>
#include <utility>

namespace CxxBench {

// typedefs
using bench_case_t = void(*)(State& state);

struct Registry {
	using bench_cases_container_t = std::vector<std::pair<std::string, bench_case_t>>;
	bench_cases_container_t bench_cases;

	void add(std::string name, bench_case_t bench_case) {
		bench_cases.push_back(std::make_pair(name, bench_case));
	}

	bench_cases_container_t::iterator begin() { return bench_cases.begin(); }
	bench_cases_container_t::iterator end() { return bench_cases.end(); }
};

} // end CxxBench