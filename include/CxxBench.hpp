#ifndef CXX_BENCH_HPP
#define CXX_BENCH_HPP

#include "CxxBench/macros.hpp"
#include "CxxBench/state.hpp"
#include "CxxBench/registry.hpp"

namespace CxxBench {
	// global variables
	extern Registry registry; // stores all benchmark cases
};

#include "CxxBench/auto_registrar.hpp"
#include "CxxBench/Timer/macros.hpp"

#endif