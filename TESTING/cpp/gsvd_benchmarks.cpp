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

#include "lapack.hpp"
#include "tools.hpp"
#include "xGGSVD3.hpp"
#include "xGGQRCS.hpp"

#include <algorithm>
#include <boost/assert.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <benchmark/benchmark.h>
#include <cassert>
#include <cstdint>
#include <limits>
#include <random>
#include <utility>



namespace ublas = boost::numeric::ublas;
namespace tools = lapack::tools;
namespace ggqrcs = lapack::ggqrcs;
namespace ggsvd3 = lapack::ggsvd3;

using Integer = lapack::integer_t;



template<template<typename Number> class Solver>
void gsvd_benchmark(benchmark::State& state) {
	using Number = float;
	using Real = float;

	auto real_nan = tools::not_a_number<Real>::value;
	auto dummy = Number{};
	auto m = static_cast<std::size_t>(state.range(0));
	auto n = static_cast<std::size_t>(state.range(1));
	auto p = static_cast<std::size_t>(state.range(0));
	auto r = std::min(m + p, n);
	auto gen = std::mt19937();

	gen.discard(1u << 17);

	auto theta_dist = ggqrcs::ThetaDistribution<Real>(0u);
	auto k = std::min( {m, p, r, m + p - r} );
	auto theta = ublas::vector<Real>(k, real_nan);

	std::generate(
		theta.begin(), theta.end(),
		[&gen, &theta_dist](){ return theta_dist(gen); }
	);

	auto min_log_cond_X = Real{0};
	auto max_log_cond_X = static_cast<Real>(std::numeric_limits<Real>::digits/2);
	auto log_cond_dist =
		std::uniform_real_distribution<Real>(min_log_cond_X, max_log_cond_X);
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
	auto solver = Solver<Number>(m, n, p);

	for(auto _ : state)
	{
		state.PauseTiming();
		solver.A = A;
		solver.B = B;
		state.ResumeTiming();

		auto ret = solver();

		BOOST_VERIFY(ret == 0);
	}
}

template<typename T> using xGGQRCS = ggqrcs::Caller<T>;
template<typename T> using xGGSVD3 = ggsvd3::Caller<T>;

BENCHMARK_TEMPLATE(
	gsvd_benchmark, xGGQRCS
)->Ranges({{8, 512}, {8, 512}});

BENCHMARK_TEMPLATE(
	gsvd_benchmark, xGGSVD3
)->Ranges({{8, 512}, {8, 512}});

BENCHMARK_MAIN();
