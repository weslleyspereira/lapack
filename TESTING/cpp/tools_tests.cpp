/*
 * Copyright (c) 2016, 2019, 2020 Christoph Conrads
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

#include <algorithm>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/test/unit_test.hpp>
#include <limits>
#include <cmath>
#include <random>


namespace ublas = boost::numeric::ublas;

using types = lapack::supported_types;
using namespace lapack::tools;



BOOST_AUTO_TEST_CASE_TEMPLATE(measure_isometry_test_simple, Number, types)
{
	for(auto m = std::size_t{0}; m < 5; ++m)
	{
		for(auto n = std::size_t{0}; n <= m; ++n)
		{
			auto A = ublas::matrix<Number, ublas::column_major>(m, n, 0);

			for(auto i = std::size_t{0}; i < n; ++i)
			{
				A(i,i) = std::pow(Number{-1}, Number(i));
			}

			BOOST_CHECK_EQUAL( 0, measure_isometry(A) );
		}
	}
}


BOOST_AUTO_TEST_CASE_TEMPLATE(measure_isometry_test_4by2, Number, types)
{
	auto A = ublas::matrix<Number, ublas::column_major>(4, 2);

	// column-major order is required for this statement to work
	std::iota( A.data().begin(), A.data().end(), 1u );

	auto I = ublas::identity_matrix<Number>(2);
	auto AT_A = ublas::matrix<Number>(2, 2);

	AT_A(0,0) = 30; AT_A(0,1) = 70;
	AT_A(1,0) = 70; AT_A(1,1) =174;

	auto expected_result = ublas::norm_frobenius(AT_A - I);

	BOOST_CHECK_EQUAL( expected_result, measure_isometry(A) );
}


BOOST_AUTO_TEST_CASE_TEMPLATE(xGEMM_test, Number, types)
{
	using Real = typename real_from<Number>::type;
	using Matrix = ublas::matrix<Number, ublas::column_major>;

	auto gen = std::minstd_rand();

	for(auto m = std::size_t{2}; m < 100; m += 10)
	{
		for(auto n = std::size_t{1}; n <= 2*m; n += 10)
		{
			auto k = std::size_t{m/2};
			auto cond = Real{1e3};
			auto A = make_matrix_like(Number(), m, k, cond, &gen);
			auto B = make_matrix_like(Number(), k, n, cond, &gen);
			auto C = ublas::prod(A, B);
			auto D = Matrix(m, n);
			auto alpha = Number{1};
			auto beta = Number{0};

			lapack::xGEMM(
				'N', 'N', m, n, k,
				alpha, &A(0,0), m, &B(0,0), k, beta, &D(0,0), m
			);

			auto eps = std::numeric_limits<Real>::epsilon();
			auto norm_C = ublas::norm_frobenius(C);
			BOOST_CHECK_LE( ublas::norm_frobenius(C-D), m*n*norm_C * eps );
		}
	}
}
