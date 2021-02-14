/*
 * Copyright (c) 2020 Christoph Conrads
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of the copyright holders nor the
 *    names of its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include "config.hpp"
#include "lapack.hpp"
#include "tools.hpp"
#include "xGGSVD3.hpp"
#include "xGGQRCS.hpp"

#include <algorithm>
#include <boost/assert.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <limits>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>



namespace ublas = boost::numeric::ublas;
namespace tools = lapack::tools;
namespace ggqrcs = lapack::ggqrcs;
namespace ggsvd3 = lapack::ggsvd3;


struct CpuClock
{
	using duration = std::chrono::duration<std::intmax_t, std::nano>;
	using time_point = std::chrono::time_point<CpuClock>;

	static time_point now()
	{
		auto tp = timespec{-1, -1};
		auto ret = clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tp);

		BOOST_VERIFY(ret == 0);

		auto second2nanosecond = std::intmax_t{1000000000};
		auto t =
			second2nanosecond * std::intmax_t{tp.tv_sec}
			+ std::intmax_t{tp.tv_nsec}
		;

		return time_point(duration(t));
	}
};



using Integer = lapack::integer_t;
using Clock = CpuClock;
using Duration = std::chrono::duration<double, std::milli>;



template<
	typename Number,
	class Prng,
	class Matrix = ublas::matrix<Number, ublas::column_major>
>
std::pair<Matrix, Matrix> make_matrix_pair(
	std::size_t m, std::size_t n, std::size_t p, Prng* p_gen
)
{
	using Real = typename tools::real_from<Number>::type;

	auto real_nan = tools::not_a_number<Real>::value;
	auto dummy = Number{};
	auto& gen = *p_gen;
	// maximize the number of singular values to compute
	auto r = std::min(m, p);
	auto k = std::min( {m, p, r, m + p - r} );

	auto option_dist = std::uniform_int_distribution<unsigned>(0, 2);
	auto min_log_cond_X = Real{0};
	auto max_log_cond_X = static_cast<Real>(std::numeric_limits<Real>::digits/2);
	auto log_cond_dist =
		std::uniform_real_distribution<Real>(min_log_cond_X, max_log_cond_X);

	auto option = option_dist(gen);
	auto theta_dist = ggqrcs::ThetaDistribution<Real>(option);
	auto theta = ublas::vector<Real>(k, real_nan);

	std::generate(
		theta.begin(), theta.end(),
		[&gen, &theta_dist](){ return theta_dist(gen); }
	);

	auto log_cond_X = log_cond_dist(gen);
	auto cond_X = std::pow(Real{2}, log_cond_X);
	auto X = tools::make_matrix_like(dummy, r, n, cond_X, &gen);
	auto U1 = tools::make_isometric_matrix_like(dummy, m, m, &gen);
	auto U2 = tools::make_isometric_matrix_like(dummy, p, p, &gen);
	auto ds = ggqrcs::assemble_diagonals_like(dummy, m, p, r, theta);
	auto D1 = ds.first;
	auto D2 = ds.second;
	auto A = ggqrcs::assemble_matrix(U1, D1, X);
	auto B = ggqrcs::assemble_matrix(U2, D2, X);

	return std::make_pair(A, B);
}



template<
	typename Number,
	template<typename T> class Solver,
	class Prng
>
std::size_t compute_num_benchmark_runs(
	std::size_t m, std::size_t n, std::size_t p,
	bool compute_matrices_p,
	Prng* p_gen
)
{
	auto t_min = std::chrono::seconds(1);
	auto t = Duration(0);
	auto ab = make_matrix_pair<Number>(m, n, p, p_gen);
	auto num_runs = std::size_t{0};

	for(auto i = std::size_t{0}; i < 32 && t < t_min; ++i)
	{
		auto cm_p = compute_matrices_p;
		auto solver = Solver<Number>(m, n, p, cm_p, cm_p, cm_p);

		for(auto j = std::size_t{0}; j < std::size_t{1}<<i; ++j)
		{
			solver.A = ab.first;
			solver.B = ab.second;
			auto t_0 = Clock::now();
			auto ret = solver();
			auto t_1 = Clock::now();

			t += t_1 - t_0;

			BOOST_VERIFY(ret == 0);
		}

		num_runs += std::size_t{1} << i;
	}

	return num_runs;
}



template<
	typename Number,
	template<typename T> class Solver
>
struct workspace_size
{
	static std::size_t get(const Solver<Number>& solver)
	{
		return solver.work.size();
	}
};

template<
	typename Real,
	template<typename T> class Solver
>
struct workspace_size<std::complex<Real>, Solver>
{
	static std::size_t get(const Solver<std::complex<Real>>& solver)
	{
		return solver.work.size() + solver.rwork.size();
	}
};

template<
	typename Number,
	template<typename T> class Solver
>
std::size_t get_workspace_size(const Solver<Number>& solver)
{
	return workspace_size<Number, Solver>::get(solver);
}


/**
 * @return Optimal workspace size, accumulated run-time.
 */
template<
	typename Number,
	template<typename T> class Solver,
	class Matrix = ublas::matrix<Number, ublas::column_major>
>
std::pair<Duration, std::size_t> run_gsvd(
	const std::vector<Matrix>& as,
	const std::vector<Matrix>& bs,
	bool compute_matrices_p
)
{
	BOOST_VERIFY(as.size() == bs.size());
	BOOST_VERIFY(as.size() > 0);

	auto m = as[0].size1();
	auto n = as[0].size2();
	auto p = bs[0].size1();
	auto cm_p = compute_matrices_p;
	auto solver = Solver<Number>(m, n, p, cm_p, cm_p, cm_p);
	auto t = Duration(0);

	for(auto i = std::size_t{0}; i < as.size(); ++i)
	{
		solver.A = as[i];
		solver.B = bs[i];
		auto t_0 = Clock::now();
		auto ret = solver();
		auto t_1 = Clock::now();

		BOOST_VERIFY(ret == 0);

		t += t_1 - t_0;
	}

	return std::make_pair(t, get_workspace_size(solver));
}


template<
	typename Number,
	template<typename T> class Solver
>
std::tuple<std::size_t, Duration, std::size_t> benchmark_gsvd(
	std::size_t m, std::size_t n, std::size_t p, unsigned seed,
	bool compute_matrices_p
)
{
	using Matrix = ublas::matrix<Number, ublas::column_major>;

	auto gen = std::mt19937(seed);

	gen.discard(1ul << 17);

	auto num_iterations_guess = compute_num_benchmark_runs<Number, Solver>(
		m, n, p, compute_matrices_p, &gen);
	auto num_iterations = std::max(
		num_iterations_guess+1, std::size_t{100}
	);
	auto as = std::vector<Matrix>(num_iterations);
	auto bs = std::vector<Matrix>(num_iterations);

	for(auto i = std::size_t{0}; i < num_iterations; ++i)
	{
		auto ab = make_matrix_pair<Number>(m, n, p, &gen);
		as[i] = ab.first;
		bs[i] = ab.second;
	}

	auto ret = run_gsvd<Number, Solver>(as, bs, compute_matrices_p);

	return std::make_tuple(num_iterations, ret.first, ret.second);
}



template<
	typename Number,
	template<typename T> class Solver,
	bool solver_built_p
>
struct BenchmarkGsvd
{
	static void run(
		const std::string& solver,
		std::size_t m, std::size_t n, std::size_t p, unsigned seed,
		bool compute_matrices_p)
	{
		auto idw =
			benchmark_gsvd<Number, Solver>(m, n, p, seed, compute_matrices_p);
		auto num_iterations = std::get<0>(idw);
		auto duration = std::get<1>(idw);
		auto workspace_size_bytes = std::get<2>(idw);
		auto t_accumulated = duration.count();
		auto t_per_sample = duration.count() / num_iterations;

		std::printf(
			"%8s  %9s  %3zu %3zu %3zu  %9zu  %8.2e  %6zu %8.2e\n",
			solver.c_str(),
			compute_matrices_p ? "Y" : "N",
			m, n, p,
			workspace_size_bytes,
			t_per_sample, num_iterations, t_accumulated
		);
	}
};

template<
	typename Number,
	template<typename T> class Solver>
struct BenchmarkGsvd<Number, Solver, false>
{
	static void run(
		const std::string&,
		std::size_t, std::size_t, std::size_t, unsigned, bool)
	{
	}
};



int main()
{
	// print column headings
	std::printf(
		"%8s  %9s  %3s %3s %3s  %9s  %8s  %6s %8s\n",
		"Solver", "Matrices?", "m", "n", "p",
		"Workspace",
		"t_CPU/ms", "Iters", "t_CPU_total/ms"
	);

	for(auto m = std::size_t{8}; m <= 128; m += (m == 8) ? 8 : 16)
	{
		for(auto n = std::size_t{8}; n <= 128; n += (n == 8) ? 8 : 16)
		{
			auto p = m;
			auto seed = 1u;

			BenchmarkGsvd<
				float, ggqrcs::Caller, lapack::BUILD_SINGLE_P
			>::run("SGGQRCS", m, n, p, seed, true);
			BenchmarkGsvd<
				float, ggsvd3::Caller, lapack::BUILD_SINGLE_P
			>::run("SGGSVD3", m, n, p, seed, true);

			BenchmarkGsvd<
				double, ggqrcs::Caller, lapack::BUILD_DOUBLE_P
			>::run("DGGQRCS", m, n, p, seed, true);
			BenchmarkGsvd<
				double, ggsvd3::Caller, lapack::BUILD_DOUBLE_P
			>::run("DGGSVD3", m, n, p, seed, true);

			BenchmarkGsvd<
				std::complex<float>, ggqrcs::Caller, lapack::BUILD_COMPLEX_P
			>::run("CGGQRCS", m, n, p, seed, true);
			BenchmarkGsvd<
				std::complex<float>, ggsvd3::Caller, lapack::BUILD_COMPLEX_P
			>::run("CGGSVD3", m, n, p, seed, true);

			BenchmarkGsvd<
				std::complex<double>, ggqrcs::Caller, lapack::BUILD_COMPLEX16_P
			>::run("ZGGQRCS", m, n, p, seed, true);
			BenchmarkGsvd<
				std::complex<double>, ggsvd3::Caller, lapack::BUILD_COMPLEX16_P
			>::run("ZGGSVD3", m, n, p, seed, true);


			BenchmarkGsvd<
				float, ggqrcs::Caller, lapack::BUILD_SINGLE_P
			>::run("SGGQRCS", m, n, p, seed, false);
			BenchmarkGsvd<
				float, ggsvd3::Caller, lapack::BUILD_SINGLE_P
			>::run("SGGSVD3", m, n, p, seed, false);

			BenchmarkGsvd<
				double, ggqrcs::Caller, lapack::BUILD_DOUBLE_P
			>::run("DGGQRCS", m, n, p, seed, false);
			BenchmarkGsvd<
				double, ggsvd3::Caller, lapack::BUILD_DOUBLE_P
			>::run("DGGSVD3", m, n, p, seed, false);

			BenchmarkGsvd<
				std::complex<float>, ggqrcs::Caller, lapack::BUILD_COMPLEX_P
			>::run("CGGQRCS", m, n, p, seed, false);
			BenchmarkGsvd<
				std::complex<float>, ggsvd3::Caller, lapack::BUILD_COMPLEX_P
			>::run("CGGSVD3", m, n, p, seed, false);

			BenchmarkGsvd<
				std::complex<double>, ggqrcs::Caller, lapack::BUILD_COMPLEX16_P
			>::run("ZGGQRCS", m, n, p, seed, false);
			BenchmarkGsvd<
				std::complex<double>, ggsvd3::Caller, lapack::BUILD_COMPLEX16_P
			>::run("ZGGSVD3", m, n, p, seed, false);
		}
	}
}
