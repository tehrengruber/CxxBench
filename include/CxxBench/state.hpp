// boost accumulators
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

// boost property tree
#include <boost/property_tree/ptree.hpp>

#include <chrono>
#include <unordered_map>

#include "Timer/rdtsc_timer.hpp"

using boost::property_tree::ptree;
using namespace boost::accumulators;

namespace CxxBench {

struct State;
#include <ostream>
std::ostream& operator<<(std::ostream& stream, const CxxBench::State& state);

struct State {
	// timer
	util::rdtsc_timer timer;

	// attributes
	std::unordered_map<std::string, std::string> attrs;

	size_t iterations=0;

	// accumulators
	accumulator_set<long double, stats<tag::sum, tag::min, tag::max, tag::mean, tag::variance(lazy)> > time_accumulator;
	accumulator_set<long double, stats<tag::sum, tag::min, tag::max, tag::mean> > cycle_accumulator;

	// return true until a predefined number of seconds has elapsed
	bool keep_running(bool output=true) {
		if (output) {
			std::cout << *this;
		}
		return sum(time_accumulator) < 5 && iterations < 100000u;
	}

	void start() {
		iterations++;
	}

	void stop() {
		time_accumulator(timer.sec());
		cycle_accumulator(timer.cycles());
	}

	template <typename T>
	void set_attr(const std::string& val, const T& attr) {
		std::stringstream ss;
		ss << attr;
		attrs.erase(val);
		attrs.emplace(val, ss.str());
	}

	ptree properties() {
		ptree pt;

		pt.put("iterations", iterations);

		// time
		ptree time_pt;
		time_pt.put<double>("min", min(time_accumulator));
		time_pt.put<double>("max", max(time_accumulator));
		time_pt.put<double>("mean", mean(time_accumulator));
		time_pt.put<double>("variance", variance(time_accumulator));
		time_pt.put<double>("stddev", std::sqrt(variance(time_accumulator)));
		pt.put_child("time", time_pt);

		// cycles
		ptree cycles_pt;
		cycles_pt.put<double>("min", min(cycle_accumulator));
		cycles_pt.put<double>("max", max(cycle_accumulator));
		cycles_pt.put<double>("mean", mean(cycle_accumulator));
		pt.put_child("cycles", cycles_pt);

		// attributes
		ptree attrs_pt;
		for (auto& attr : attrs) {
			attrs_pt.put(attr.first, attr.second);
		}

		if (attrs.size() > 0)
			pt.put_child("attrs", attrs_pt);

		return pt;
	}

	/*static std::string header() {
		std::stringstream header;

		auto additional_attr_headers = argument_t::header_attributes();
		for (std::string& attr : additional_attr_headers) {
			header << std::setw(15) << attr << " ";
		}

		header << std::setw(12) << "iterations" << " "
			   << std::setw(12) << "mean" << " "
			   << std::setw(12) << "min"  << " "
			   << std::setw(12) << "max";
		return header.str();
	}*/
};

std::ostream& operator<<(std::ostream& stream, const CxxBench::State& state) {
	// no output in the first iteration
	if (state.iterations==0)
		return stream;

	// print header after first iteration
	if (state.iterations==1) {
		for (auto& row : state.attrs) {
			std::cout << std::setw(15) << row.first << " ";
		}

		std::cout << std::setw(12) << "iterations" << " "
				  << std::setw(10) << "mean" << " "
				  << std::setw(10) << "min"  << " "
				  << std::setw(10) << "max" << std::endl;
	}

	printf("\r");

	for (auto& row : state.attrs) {
		stream << std::setw(15) << row.second << " ";
	}

	stream << std::setw(12) << state.iterations << " ";

	if (mean(state.time_accumulator) < 0.001)
		stream << std::scientific;
	else
		stream << std::fixed;

	stream
		<< std::setw(12) << mean(state.time_accumulator) << " "
		<< std::setw(12) << min(state.time_accumulator) << " "
		<< std::setw(12) << max(state.time_accumulator)
		/*<< std::defaultfloat*/;

	stream.unsetf(std::ios_base::floatfield); // fix for gcc-4.9 instead of using std::defaultfloat

	if (mean(state.time_accumulator) < 0.001)
		stream << std::scientific;

	//stream << std::endl;
	stream << std::flush;

	return stream;
}

} // end CxxBench