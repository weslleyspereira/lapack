/*
 * Copyright (c) 2020 Christoph Conrads (https://christoph-conrads.name)
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
#include "tools.hpp"

#include <cstddef>
#include <boost/assert.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/test/unit_test.hpp>

using types = lapack::supported_types;
using Integer = lapack::integer_t;

namespace ublas = boost::numeric::ublas;
namespace tools = lapack::tools;


//
// certain xORCSD2BY1, xUNCSD2BY1 problems must have been fixed for xGGQRCS to
// work; test for these problems here.
//

template<
	typename Number,
	typename std::enable_if<
		std::is_fundamental<Number>::value, int
	>::type* = nullptr
>
void xUNCSD2BY1_regression_20200420_impl(Number)
{
	using Real = typename tools::real_from<Number>::type;
	using Matrix = ublas::matrix<Number, ublas::column_major>;

	auto m = std::size_t{130}; // 160
	auto n = std::size_t{131}; // 179
	auto p = std::size_t{130}; // 150
	auto ldq = m + p;
	auto gen = std::minstd_rand();
	auto Q = tools::make_isometric_matrix_like(Real{0}, ldq, n, &gen);

	BOOST_VERIFY( tools::is_almost_isometric(Q) );

	auto nan = tools::not_a_number<Number>::value;
	auto real_nan = tools::not_a_number<Real>::value;
	auto theta = ublas::vector<Real>(n, real_nan);
	auto U1 = Matrix(m, m, nan);
	auto U2 = Matrix(p, p, nan);
	auto Qt = Matrix(n, n, nan);
	auto lwork = 32 * ldq;
	auto work = ublas::vector<Number>(lwork, nan);
	auto iwork = ublas::vector<Integer>(m+p, -1);
	auto ret = lapack::xUNCSD2BY1(
		'Y', 'Y', 'N',
		m + p, m, n,
		&Q(0, 0), ldq, &Q(m, 0), ldq,
		&theta(0),
		&U1(0,0), m, &U2(0,0), p, &Qt(0,0), n,
		&work(0), lwork, &iwork(0)
	);

	BOOST_REQUIRE_EQUAL( ret, 0 );
	BOOST_CHECK( tools::is_almost_isometric(U1) );
	BOOST_CHECK( tools::is_almost_isometric(U2) );
}

template<
	typename Number,
	typename std::enable_if<
		!std::is_fundamental<Number>::value, int
	>::type* = nullptr
>
void xUNCSD2BY1_regression_20200420_impl(Number)
{
	using Real = typename tools::real_from<Number>::type;
	using Matrix = ublas::matrix<Number, ublas::column_major>;

	auto m = std::size_t{130}; // 160
	auto n = std::size_t{131}; // 179
	auto p = std::size_t{130}; // 150
	auto ldq = m + p;
	auto gen = std::minstd_rand();
	auto Q = tools::make_isometric_matrix_like(Number{0}, ldq, n, &gen);

	BOOST_VERIFY( tools::is_almost_isometric(Q) );

	auto nan = tools::not_a_number<Number>::value;
	auto real_nan = tools::not_a_number<Real>::value;
	auto theta = ublas::vector<Real>(n, real_nan);
	auto U1 = Matrix(m, m, nan);
	auto U2 = Matrix(p, p, nan);
	auto Qt = Matrix(n, n, nan);
	auto lwork = 32 * ldq;
	auto work = ublas::vector<Number>(lwork, nan);
	auto lrwork = 32 * ldq;
	auto rwork = ublas::vector<Real>(lrwork, real_nan);
	auto iwork = ublas::vector<Integer>(m+p, -1);
	auto ret = lapack::xUNCSD2BY1(
		'Y', 'Y', 'N',
		m + p, m, n,
		&Q(0, 0), ldq, &Q(m, 0), ldq,
		&theta(0),
		&U1(0,0), m, &U2(0,0), p, &Qt(0,0), n,
		&work(0), lwork, &rwork(0), lrwork, &iwork(0)
	);

	BOOST_REQUIRE_EQUAL( ret, 0 );
	BOOST_CHECK( tools::is_almost_isometric(U1) );
	BOOST_CHECK( tools::is_almost_isometric(U2) );
}

BOOST_AUTO_TEST_CASE_TEMPLATE(xUNCSD2BY1_regression_20200420, Number, types)
{
	xUNCSD2BY1_regression_20200420_impl(Number{});
}



template<
	typename Number,
	typename std::enable_if<
		std::is_fundamental<Number>::value, int
	>::type* = nullptr
>
void xUNCSD2BY1_test_workspace_query_with_lrwork_impl(Number)
{
	// there is no lrwork parameter for real-valued 2-by-1 CSD
	BOOST_CHECK_EQUAL( true, true ); // silence Boost warnings
}

template<
	typename Number,
	typename std::enable_if<
		!std::is_fundamental<Number>::value, int
	>::type* = nullptr
>
void xUNCSD2BY1_test_workspace_query_with_lrwork_impl(Number)
{
	using Real = typename tools::real_from<Number>::type;

	// The workspace sizes are returned as floating-point values. Thus,
	// the dimensions should not be too large if you want to check the
	// workspace query result.
	auto m = std::size_t{(1u<<8)+1};
	auto n = std::size_t{(1u<<16)+1};
	auto p = std::size_t{(1u<<19)-1};
	auto nan = tools::not_a_number<Number>::value;
	auto real_nan = tools::not_a_number<Real>::value;
	auto dummy = nan;
	auto real_dummy = real_nan;
	auto ldq = 2 * (m + p) + 1;
	auto ldu1 = 4 * m + 7;
	auto ldu2 = 3 * p + 5;
	auto ldqt = 2 * n + 13;
	auto work = nan;
	auto lwork = 1;
	auto rwork = real_nan;
	auto lrwork = -1;
	auto iwork = Integer{-128};
	auto ret = lapack::xUNCSD2BY1(
		'Y', 'Y', 'Y',
		m + p, m, n,
		&dummy, ldq, &dummy, ldq,
		&real_dummy,
		&dummy, ldu1, &dummy, ldu2, &dummy, ldqt,
		&work, lwork, &rwork, lrwork, &iwork
	);

	auto k = std::min({ m, n, p, m+p-n });

	BOOST_REQUIRE_EQUAL( ret, 0 );
	BOOST_CHECK_GE( std::real(work), 8*k );
	BOOST_CHECK_GE( rwork, 8*k );
}

BOOST_AUTO_TEST_CASE_TEMPLATE(
	xUNCSD2BY1_test_workspace_query_with_lrwork, Number, types)
{
	xUNCSD2BY1_test_workspace_query_with_lrwork_impl(Number{});
}
