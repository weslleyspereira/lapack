/*
 * Copyright (c) 2020, Christoph Conrads (https://christoph-conrads.name)
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

#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <type_traits>

namespace tools = lapack::tools;

using Integer = lapack::integer_t;
using types = lapack::supported_types;

//
// ATTENTION
//
// One-based indexing must be used for all indices passed to xLASRTI.
//

template<
	typename Number,
	typename std::enable_if<
		std::is_fundamental<Number>::value, int
	>::type* = nullptr
>
void xLASRTI_test_simple_impl(Number)
{
	using Real = typename tools::real_from<Number>::type;

	auto xs = std::vector<Real>( { 3, 2, -1 } );
	auto is = std::vector<Integer>( { 1, 2, 3 } );
	auto ret = lapack::xLASRTI('D', 3, xs.data(), is.data());

	BOOST_REQUIRE_EQUAL( ret, 0 );

	BOOST_CHECK_EQUAL( xs[0], 3 );
	BOOST_CHECK_EQUAL( xs[1], 2 );
	BOOST_CHECK_EQUAL( xs[2], -1 );

	BOOST_CHECK_EQUAL( is[0], 1 );
	BOOST_CHECK_EQUAL( is[1], 2 );
	BOOST_CHECK_EQUAL( is[2], 3 );

	ret = lapack::xLASRTI('I', 3, xs.data(), is.data());

	BOOST_REQUIRE_EQUAL( ret, 0 );

	BOOST_CHECK_EQUAL( xs[0], 3 );
	BOOST_CHECK_EQUAL( xs[1], 2 );
	BOOST_CHECK_EQUAL( xs[2], -1 );

	BOOST_CHECK_EQUAL( is[0], 3 );
	BOOST_CHECK_EQUAL( is[1], 2 );
	BOOST_CHECK_EQUAL( is[2], 1 );
}

template<
	typename Number,
	typename std::enable_if<
		!std::is_fundamental<Number>::value, int
	>::type* = nullptr
>
void xLASRTI_test_simple_impl(Number)
{
	// complex values cannot be sorted
	BOOST_CHECK(true); // suppress Boost warnings
}



BOOST_AUTO_TEST_CASE_TEMPLATE(xLASRTI_test_simple, Number, types)
{
	xLASRTI_test_simple_impl(Number{});
}


template<
	typename Number,
	typename std::enable_if<
		std::is_fundamental<Number>::value, int
	>::type* = nullptr
>
void xLASRTI_test_shuffle_impl(Number)
{
	using Real = typename tools::real_from<Number>::type;

	// test at least up to 21 because for n<=20 insertion sort is used
	for(auto n = std::size_t{0}; n <= 30; ++n)
	{
		auto gen = std::mt19937(n+1);
		auto xs = std::vector<Real>(n);
		auto is = std::vector<Integer>(n);

		std::iota(xs.begin(), xs.end(), 0);
		std::iota(is.begin(), is.end(), 1);
		std::shuffle(is.begin(), is.end(), gen);

		// sort (increasing order)
		auto ret = lapack::xLASRTI('I', n, xs.data(), is.data());

		BOOST_REQUIRE_EQUAL( ret, 0 );

		for(auto k = std::size_t{0}; k < n; ++k)
		{
			BOOST_CHECK_EQUAL( xs[k], k );
			BOOST_CHECK_EQUAL( is[k], k+1 );
		}


		// call sort on sorted data
		ret = lapack::xLASRTI('I', n, xs.data(), is.data());

		BOOST_REQUIRE_EQUAL( ret, 0 );

		for(auto k = std::size_t{0}; k < n; ++k)
		{
			BOOST_CHECK_EQUAL( xs[k], k );
			BOOST_CHECK_EQUAL( is[k], k+1 );
		}


		// sort (decreasing order)
		ret = lapack::xLASRTI('D', n, xs.data(), is.data());

		BOOST_REQUIRE_EQUAL( ret, 0 );

		for(auto k = std::size_t{0}; k < n; ++k)
		{
			BOOST_CHECK_EQUAL( xs[k], k );
			BOOST_CHECK_EQUAL( is[k], n-k );
		}


		// sort again after underlying array changed
		std::shuffle(xs.begin(), xs.end(), gen);

		ret = lapack::xLASRTI('D', n, xs.data(), is.data());

		BOOST_REQUIRE_EQUAL( ret, 0 );

		for(auto k = std::size_t{1}; k < n; ++k)
		{
			BOOST_CHECK_GE( xs[is[k-1]-1], xs[is[k]-1] );
		}

		std::sort(is.begin(), is.end());

		// check no indices were erased
		for(auto k = std::size_t{0}; k < n; ++k)
		{
			BOOST_CHECK_EQUAL( is[k], k+1 );
		}
	}
}

template<
	typename Number,
	typename std::enable_if<
		!std::is_fundamental<Number>::value, int
	>::type* = nullptr
>
void xLASRTI_test_shuffle_impl(Number)
{
	// complex values cannot be sorted
	BOOST_CHECK(true); // suppress Boost warnings
}

BOOST_AUTO_TEST_CASE_TEMPLATE(xLASRTI_test_shuffle, Number, types)
{
	xLASRTI_test_shuffle_impl(Number{});
}


