#include <CxxBench.hpp>
#include <cmath>
#include <cstdint>
#include <x86intrin.h>
#include <cassert>

constexpr double G = 6.674e-11; 

/**
 * \param N number of bodies
 * \param x x coords
 * \param y y coords
 * \param z z coords
 * \param fx computed force f_x
 * \param fy computed force f_y 
 * \param fz computed force f_z
 */
void nbody(uint32_t N, float* x, float* y, float* z, float* fx, float* fy, float* fz) {
	for (unsigned i=0; i<N; ++i) {
		fx[i] = 0;
		fy[i] = 0;
		fz[i] = 0;
		for (unsigned j=0; j<N; ++j) {
			if (i==j)
				continue;

			float xij = x[i]-x[j];
			float yij = y[i]-y[j];
			float zij = z[i]-z[j];
	 		float rij3 = ((xij*xij)+(yij*yij))+(zij*zij);
	 		rij3 *= std::sqrt(rij3);
	 		fx[i] += xij / rij3;
	 		fy[i] += yij / rij3;
	 		fz[i] += zij / rij3;
	 	}
	 	fx[i] *= G;
		fy[i] *= G;
		fz[i] *= G;
	}
}

void nbody_vec(uint32_t N, float* x, float* y, float* z, 
				float* fx, float* fy, float* fz) {
	assert(N%8 == 0);

	// define intervals 

	for (uint32_t i=0; i<N; i++) {
		// set forces to zero
		fx[i] = 0;
		fy[i] = 0;
		fz[i] = 0;

		// load current values
		__m256 vxi = _mm256_set1_ps(*(x+i));
		__m256 vyi = _mm256_set1_ps(*(y+i));
		__m256 vzi = _mm256_set1_ps(*(z+i));

		// initialize forces with zero
		__m256 vfxi = _mm256_setzero_ps();
		__m256 vfyi = _mm256_setzero_ps();
		__m256 vfzi = _mm256_setzero_ps();

		for (uint32_t j=0; j<N; j+=8) {
			// serial fallback if i \in [j, j+7]
			if (i >= j && i < j+8) {
				for (unsigned k=j; k<j+8; ++k) {
					if (i==k)
						continue;

					float xij = x[i]-x[k];
					float yij = y[i]-y[k];
					float zij = z[i]-z[k];
			 		float rij3 = xij*xij+yij*yij+zij*zij;
			 		rij3 *= std::sqrt(rij3);
			 		fx[i] += xij / rij3;
			 		fy[i] += yij / rij3;
			 		fz[i] += zij / rij3;
			 	}
				continue;
			}

			// compute x_ij, y_ij, z_ij
			__m256 vxij = _mm256_sub_ps(vxi, _mm256_load_ps(x+j));
			__m256 vyij = _mm256_sub_ps(vyi, _mm256_load_ps(y+j));
			__m256 vzij = _mm256_sub_ps(vzi, _mm256_load_ps(z+j));

			// compute rij^3
			//  first compute rij2
			__m256 vrij3 = _mm256_add_ps(_mm256_mul_ps(vxij, vxij), _mm256_mul_ps(vyij, vyij));
			vrij3 = _mm256_add_ps(vrij3, _mm256_mul_ps(vzij, vzij));
			//  multiply with rij
			vrij3 = _mm256_mul_ps(vrij3, _mm256_sqrt_ps(vrij3));

			// compute force
			__m256 dfx = _mm256_div_ps(vxij, vrij3);
			__m256 dfy = _mm256_div_ps(vyij, vrij3);
			__m256 dfz = _mm256_div_ps(vzij, vrij3);
			for (uint32_t j=0; j<8; j++) {
				fx[i] += ((float*) &dfx)[j];
				fy[i] += ((float*) &dfy)[j];
				fz[i] += ((float*) &dfz)[j];
			}

			// compute force
			//vfxi = _mm256_add_ps(vfxi, _mm256_div_ps(vxij, vrij3));
			//vfyi = _mm256_add_ps(vfyi, _mm256_div_ps(vyij, vrij3));
			//vfzi = _mm256_add_ps(vfzi, _mm256_div_ps(vzij, vrij3));
		}

		// horizontal add forces in a naive way
		for (uint32_t j=0; j<8; j++) {
			fx[i] += ((float*) &vfxi)[j];
			fy[i] += ((float*) &vfyi)[j];
			fz[i] += ((float*) &vfzi)[j];
		}

		// multiply by G
		fx[i] *= G;
		fy[i] *= G;
		fz[i] *= G;
	}
}

namespace CxxBench {
	// global variables
	Registry registry; // stores all benchmark cases
};

int main () {
	for (auto& bench_case : CxxBench::registry) {
		CxxBench::State state;

		// print header
		std::string name = bench_case.first;
		CxxBench::bench_case_t bench_func = bench_case.second;
		std::cout << "-----------------------------------------------" << std::endl
			 	  << " Benchmark " << name << std::endl
			 	  << "-----------------------------------------------" << std::endl;


		// run benchmark
		bench_func(state);
	}
}

CXX_BENCHMARK_CASE(nbody_problem) {
	int N = 4096;
	state.set_attr("N", N);
	float* x = (float*)_mm_malloc(N*sizeof(float), 32);
	float* y = (float*)_mm_malloc(N*sizeof(float), 32);
	float* z = (float*)_mm_malloc(N*sizeof(float), 32);
	float* fx = (float*)_mm_malloc(N*sizeof(float), 32);
	float* fy = (float*)_mm_malloc(N*sizeof(float), 32);
	float* fz = (float*)_mm_malloc(N*sizeof(float), 32);

	// create particles on a regular lattice
	int nx, ny, nz;
	nx = ny = nz = std::cbrt(N);
	assert(nx*ny*nz == N);
	double hx, hy, hz;
	hx = hy = hz = 0.01/(nx-1);

	int n=0;
	for (uint32_t i=0; i<nx; ++i) {
		for (uint32_t j=0; j<ny; ++j) {
			for (uint32_t k=0; k<nz; ++k) {
				x[n] = i*hx;
				y[n] = j*hy;
				z[n] = k*hz;

				++n;
			}
		}
	}

	// run benchmark
	while(state.keep_running()) {
		CXX_BENCH_TIMER_START(state)
		nbody_vec(N, x, y, z, fx, fy, fz);
		CXX_BENCH_TIMER_STOP(state)
		state.set_attr("flops", (16*N*N+3*N)*state.iterations/sum(state.time_accumulator));
	}
}

//CXX_BENCH_CASE_ARG_TYPE(nbody_problem, N) {
//	int N = 4096;
//}
//
//CXX_BENCH_CASE_ARG_GENERATOR(nbody_problem) {
//	std::vector<arg_t> args;
//	for (int i=128; i<=4096; ++i) {
//		args.emplace_back(i);
//	}
//}
//
//struct NBody : CxxBench::Benchmark::Case {
//	CXX_BENCH_CASE_ARG_LIST(N)
//
//	struct Argument {
//		int N;
//	};
//
//	void operator(CxxBench::Benchmark::State& state) {
//
//	}
//}