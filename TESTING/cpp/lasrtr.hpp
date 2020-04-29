/*
 * Copyright (c) 2020 Christoph Conrads (https://christoph-conrads.name)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *  * Redistributions of source code must retain the above copyright
 *	notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *	notice, this list of conditions and the following disclaimer in the
 *	documentation and/or other materials provided with the distribution.
 *  * Neither the name of the copyright holders nor the
 *	names of its contributors may be used to endorse or promote products
 *	derived from this software without specific prior written permission.
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
#ifndef LAPACK_TESTS_xLASRTR_HPP
#define LAPACK_TESTS_xLASRTR_HPP

#include "lapack.hpp"
#include "tools.hpp"

#include <algorithm>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/test/unit_test.hpp>
#include <limits>
#include <random>
#include <vector>


namespace ublas = boost::numeric::ublas;

using Integer = lapack::integer_t;

template<typename Number>
void xLASRTR_test_simple_impl(Number)
{
	using Real = typename real_from<Number>::type;
	using Matrix = ublas::matrix<Number, ublas::column_major>;

	auto A = Matrix(3, 1);

	A(0,0) = 1;
	A(1,0) = 2;
	A(2,0) = 3;

	auto ipvt = std::vector<Integer>(3, -1);
	auto real_nan = not_a_number<Real>::value;
	auto work = std::vector<Real>(3, real_nan);
	auto ret = lapack::lasrtr('D', 3, 1, &A(0,0), 3, ipvt.data(), work.data());

	BOOST_REQUIRE_EQUAL( ret, 0 );

	BOOST_CHECK_EQUAL( A(0,0), 3 );
	BOOST_CHECK_EQUAL( A(1,0), 2 );
	BOOST_CHECK_EQUAL( A(2,0), 1 );

	BOOST_CHECK_EQUAL( ipvt[0], 3 );
	BOOST_CHECK_EQUAL( ipvt[1], 2 );
	BOOST_CHECK_EQUAL( ipvt[2], 1 );

	BOOST_CHECK_EQUAL( work[0], 1 );
	BOOST_CHECK_EQUAL( work[1], 2 );
	BOOST_CHECK_EQUAL( work[2], 3 );

	ret = lapack::lasrtr('I', 3, 1, &A(0,0), 3, ipvt.data(), work.data());

	BOOST_REQUIRE_EQUAL( ret, 0 );

	BOOST_CHECK_EQUAL( A(0,0), 1 );
	BOOST_CHECK_EQUAL( A(1,0), 2 );
	BOOST_CHECK_EQUAL( A(2,0), 3 );

	BOOST_CHECK_EQUAL( ipvt[0], 3 );
	BOOST_CHECK_EQUAL( ipvt[1], 2 );
	BOOST_CHECK_EQUAL( ipvt[2], 1 );

	BOOST_CHECK_EQUAL( work[0], 3 );
	BOOST_CHECK_EQUAL( work[1], 2 );
	BOOST_CHECK_EQUAL( work[2], 1 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE(xLASRTR_test_simple, Number, test_types)
{
	xLASRTR_test_simple_impl(Number{});
}



template<typename Number>
void xLASRTR_test_random_impl(Number)
{
	using Real = typename real_from<Number>::type;
	using Matrix = ublas::matrix<Number, ublas::column_major>;

	// Hint for choice of dimensions:
	// xLASRT is using introsort until dimension 20
	// choice of loop counter increment: speed up test
	for(auto m = std::size_t{1}; m <= 32; m += (m<4) ? 1 : 4)
	{
		for(auto n = std::size_t{1}; n <= 25; n += (n<4) ? 1 : 3)
		{
			BOOST_TEST_CONTEXT("m=" << m) {
			BOOST_TEST_CONTEXT("n=" << n) {

			auto gen = std::minstd_rand(m + 257*n);
			auto dist = UniformDistribution<Number>();
			auto lda = m + (m*n) % 13 + 7;
			auto nan = not_a_number<Number>::value;
			auto A = Matrix(lda, n, nan);

			for(auto j = std::size_t{0}; j < n; ++j)
			{
				for(auto i = std::size_t{0}; i < m; ++i)
				{
					A(i,j) = dist(gen);
				}
			}

			auto ipvt = std::vector<Integer>(m, -1);
			auto real_nan = not_a_number<Real>::value;
			auto work = std::vector<Real>(m, real_nan);

			// sort rows, check decreasing maximum norm values
			auto ret = lapack::lasrtr(
				'D', m, n, &A(0,0), lda, ipvt.data(), work.data()
			);

			BOOST_REQUIRE_EQUAL( ret, 0 );

			for(auto i = std::size_t{1}; i < m; ++i)
			{
				auto x = ublas::row(A, i-1);
				auto y = ublas::row(A, i);

				BOOST_CHECK_GE( ublas::norm_inf(x), ublas::norm_inf(y) );
			}

			auto finite_p = [] (Real x) { return std::isfinite(x); };

			if(m > 1)
				BOOST_REQUIRE(std::all_of(work.begin(), work.end(), finite_p));

			std::sort(ipvt.begin(), ipvt.end());

			for(auto i = std::size_t{1}; i < m; ++i)
			{
				BOOST_CHECK_EQUAL( ipvt[i-1] + 1, ipvt[i] );
			}


			// shuffle rows, sort, then check matrix was restored	
			auto B = Matrix(A);

			std::shuffle(ipvt.begin(), ipvt.end(), gen);
			lapack::lapmr(true, m, n, &A(0,0), lda, ipvt.data());

			ret = lapack::lasrtr(
				'D', m, n, &A(0,0), lda, ipvt.data(), work.data()
			);
			BOOST_REQUIRE_EQUAL( ret, 0 );

			for(auto j = std::size_t{0}; j < n; ++j)
			{
				for(auto i = std::size_t{0}; i < m; ++i)
				{
					BOOST_TEST_CONTEXT("i=" << i) {
						BOOST_TEST_CONTEXT("j=" << j) {
							BOOST_CHECK_EQUAL( A(i,j), B(i,j) );
						}
					}
				}
			}


			// check nan part was not touched
			auto C = Matrix(ublas::subrange(A, m+1, lda, 0, n));
			auto nan_p = [] (Number x) {
				return std::isnan(std::real(x)) || std::isnan(std::imag(x));
			};

			BOOST_CHECK(std::all_of(C.data().begin(), C.data().end(), nan_p));
			}
			}
		}
	}
}

BOOST_AUTO_TEST_CASE_TEMPLATE(xLASRTR_test_random, Number, test_types)
{
	xLASRTR_test_random_impl(Number{});
}

#endif
