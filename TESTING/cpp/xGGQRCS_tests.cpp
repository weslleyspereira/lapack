/*
 * Copyright (c) 2020-2021 Christoph Conrads
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
#include "xGGQRCS.hpp"
#include "xGGSVD3.hpp"

#include <boost/assert.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/test/unit_test.hpp>
#include <ctime>


namespace ggqrcs = lapack::ggqrcs;
namespace ggsvd3 = lapack::ggsvd3;
namespace tools = lapack::tools;
namespace ublas = boost::numeric::ublas;

using types = lapack::supported_types;



/**
 * This LAPACK xerbla implementation prints an error message but does not
 * terminate the program thereby allowing the calling LAPACK function to return
 * to its caller.
 *
 * @param[in] caller A string WITHOUT ZERO TERMINATOR
 * @param[in] caller_len The length of the string referenced by caller
 */
extern "C" void xerbla_(
	const char* caller, int* p_info, std::size_t caller_len)
{
	BOOST_VERIFY( caller != nullptr );
	BOOST_VERIFY( p_info != nullptr );

	// "sz" prefix taken from hungarian notation (zero-terminated string)
	char szCaller[80];
	auto num_bytes_to_copy = std::min(sizeof(szCaller)-1, caller_len);

	std::memset(szCaller, 0, sizeof(szCaller));
	std::strncpy(szCaller, caller, num_bytes_to_copy);
	std::fprintf(
		stderr, "%s: parameter %d has illegal value\n", szCaller, *p_info
	);
}


BOOST_AUTO_TEST_CASE_TEMPLATE(xGGQRCS_test_simple_1x2, Number, types)
{
	auto m = std::size_t{1};
	auto n = std::size_t{2};
	auto p = std::size_t{1};
	auto caller = ggqrcs::Caller<Number>(m, n, p);
	auto A = caller.A;
	auto B = caller.B;

	A(0,0) = 1; A(0,1) = +0;
	B(0,0) = 1; B(0,1) = -1;

	caller.A = A;
	caller.B = B;

	auto ret = caller();
	check_results(ret, A, B, caller);

	BOOST_CHECK_EQUAL( caller.rank, 2 );
}


BOOST_AUTO_TEST_CASE_TEMPLATE(xGGQRCS_test_simple_2x2, Number, types)
{
	auto m = std::size_t{2};
	auto n = std::size_t{2};
	auto p = std::size_t{2};
	auto caller = ggqrcs::Caller<Number>(m, n, p);
	auto A = caller.A;
	auto B = caller.B;

	A(0,0) = 1;
	B(1,1) = 1;

	caller.A = A;
	caller.B = B;

	auto ret = caller();
	check_results(ret, A, B, caller);

	BOOST_CHECK_EQUAL( caller.rank, 2 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE(xGGQRCS_test_simple_1_2_3, Number, types)
{
	auto m = std::size_t{1};
	auto n = std::size_t{2};
	auto p = std::size_t{3};
	auto caller = ggqrcs::Caller<Number>(m, n, p);
	auto A = caller.A;
	auto B = caller.B;

	A(0,0) = 1; A(0,1) = 1;
	            B(0,1) = 100;

	caller.A = A;
	caller.B = B;

	auto ret = caller();

	check_results(ret, A, B, caller);
	BOOST_CHECK_EQUAL( caller.rank, 2 );
}



BOOST_AUTO_TEST_CASE_TEMPLATE(
	xGGQRCS_test_workspace_size_check, Number, lapack::supported_real_types
)
{
	using Real = typename tools::real_from<Number>::type;

	auto m = 17;
	auto n = 13;
	auto p = 8;
	auto rank = -1;
	auto swapped_p = false;
	auto a = Real{1};
	auto b = Real{2};
	auto alpha = Real{-1};
	auto beta = Real{-1};
	auto u1 = Real{0};
	auto u2 = Real{0};
	auto tol = Real{0};
	auto work = Real{0};
	auto iwork = -1;
	auto info = lapack::xGGQRCS(
		'Y', 'Y', 'Y', m, n, p, &rank, &swapped_p,
		&a, m, &b, p,
		&alpha, &beta,
		&u1, m, &u2, p,
		&tol,
		&work, 1, &iwork
	);

	BOOST_CHECK_EQUAL(info, -20);
}



BOOST_AUTO_TEST_CASE_TEMPLATE(
	xGGQRCS_test_workspace_size_check_complex, Number, lapack::supported_complex_types
)
{
	using Real = typename tools::real_from<Number>::type;

	auto m = 17;
	auto n = 13;
	auto p = 8;
	auto rank = -1;
	auto swapped_p = false;
	auto a = Number{1};
	auto b = Number{2};
	auto alpha = Real{-1};
	auto beta = Real{-1};
	auto tol = Real{-1};
	auto u1 = Number{0};
	auto u2 = Number{0};
	auto work = Number{0};
	auto rwork = Real{0};
	auto iwork = -1;
	auto info = lapack::xGGQRCS(
		'Y', 'Y', 'Y', m, n, p, &rank, &swapped_p,
		&a, m, &b, p,
		&alpha, &beta,
		&u1, m, &u2, p,
		&tol,
		&work, 1, &rwork, 1024, &iwork
	);

	BOOST_CHECK_EQUAL(info, -20);

	info = lapack::xGGQRCS(
		'Y', 'Y', 'Y', m, n, p, &rank, &swapped_p,
		&a, m, &b, p,
		&alpha, &beta,
		&u1, m, &u2, p,
		&tol,
		&work, 1024, &rwork, 1, &iwork
	);

	BOOST_CHECK_EQUAL(info, -22);
}




// this test does not pass with row sorting and no matrix scaling
BOOST_AUTO_TEST_CASE(xGGQRCS_test_matrix_scaling)
{
	using Number = float;

	auto m = std::size_t{1};
	auto n = std::size_t{2};
	auto p = std::size_t{10};
	auto caller = ggqrcs::Caller<Number>(m, n, p);
	auto A = caller.A;
	auto B = caller.B;

	A(0,0) = -8.519847412e+02; A(0,1) = +6.469862671e+02;
	B(0,0) = +5.485938125e+05; B(0,1) = -4.166526250e+05;
	B(1,0) = +1.846850781e+05; B(1,1) = -1.402660781e+05;
	B(2,0) = +5.322575625e+05; B(2,1) = -4.042448438e+05;
	B(3,0) = -1.630551465e+04; B(3,1) = +1.238360352e+04;
	B(4,0) = -1.286453438e+05; B(4,1) = +9.770555469e+04;
	B(5,0) = -1.323287812e+05; B(5,1) = +1.005026797e+05;
	B(6,0) = +5.681228750e+05; B(6,1) = -4.314841250e+05;
	B(7,0) = -3.107875312e+05; B(7,1) = +2.360408594e+05;
	B(8,0) = +1.456551719e+05; B(8,1) = -1.106233281e+05;
	B(9,0) = +1.365355156e+05; B(9,1) = -1.036972344e+05;

	caller.A = A;
	caller.B = B;

	auto ret = caller();

	check_results(ret, A, B, caller);
}



BOOST_TEST_DECORATOR(* boost::unit_test::expected_failures(1))
BOOST_AUTO_TEST_CASE(xGGQRCS_test_conditional_backward_stability)
{
	using Number = float;
	using Real = typename tools::real_from<Number>::type;
	using Matrix = ublas::matrix<Number, ublas::column_major>;

	constexpr auto eps = std::numeric_limits<Real>::epsilon();
	auto tol = [] (const Matrix& a) {
		return std::max(a.size1(), a.size2()) * ublas::norm_frobenius(a) * eps;
	};

	auto m = std::size_t{2};
	auto n = std::size_t{2};
	auto p = std::size_t{2};
	auto r = std::size_t{2};
	auto A = Matrix(m, n);
	auto B = Matrix(p, n);

	A(0,0) = +4.663013916e+02; A(0,1) = +4.046628418e+02;
	A(1,0) = +3.062543457e+03; A(1,1) = +2.648934082e+03;
	B(0,0) = -2.966550887e-01; B(0,1) = -2.563934922e-01;
	B(1,0) = +7.012547851e-01; B(1,1) = +6.062732935e-01;

	// compute GSVD of A, B and fail
	{
		auto caller = ggqrcs::Caller<Number>(m, n, p);

		caller.A = A;
		caller.B = B;

		auto ret = caller();
		BOOST_REQUIRE_EQUAL( ret, 0 );

		check_results(ret, A, B, caller);

		auto X = copy_X(caller);
		auto ds = ggqrcs::assemble_diagonals_like(
			Number{}, m, p, r, caller.swapped_p, caller.alpha, caller.beta
		);
		auto& D1 = ds.first;
		auto& D2 = ds.second;
		auto almost_A = ggqrcs::assemble_matrix(caller.U1, D1, X);
		auto almost_B = ggqrcs::assemble_matrix(caller.U2, D2, X);

		BOOST_CHECK_LE(ublas::norm_frobenius(A - almost_A), tol(A));
		// should fail
		BOOST_CHECK_LE(ublas::norm_frobenius(B - almost_B), tol(B));
	}

	// try again with norm(A) = norm(B)
	{
		auto w = std::ldexp(Real{1}, 12);
		auto caller = ggqrcs::Caller<Number>(m, n, p);

		caller.A = A;
		caller.B = w * B;

		auto ret = caller();
		BOOST_REQUIRE_EQUAL( ret, 0 );

		auto X = copy_X(caller);
		auto ds = ggqrcs::assemble_diagonals_like(
			Number{}, m, p, r, caller.swapped_p, caller.alpha, caller.beta
		);
		auto& D1 = ds.first;
		auto& D2 = ds.second;
		auto almost_A = ggqrcs::assemble_matrix(caller.U1, D1, X);
		auto almost_B = ggqrcs::assemble_matrix(caller.U2, D2, X);

		BOOST_CHECK_LE(ublas::norm_frobenius(A - almost_A), tol(A));
		BOOST_CHECK_LE(ublas::norm_frobenius(w*B - almost_B), tol(w*B));
	}
}



BOOST_TEST_DECORATOR(* boost::unit_test::expected_failures(1))
BOOST_AUTO_TEST_CASE(xGGQRCS_test_singular_accuracy_vs_radians_accuracy)
{
	using Number = float;
	using Real = typename tools::real_from<Number>::type;
	using Matrix = ublas::matrix<Number, ublas::column_major>;

	constexpr auto eps = std::numeric_limits<Real>::epsilon();
	auto m = std::size_t{2};
	auto n = std::size_t{1};
	auto p = std::size_t{2};
	auto r = std::size_t{1};
	auto A = Matrix(m, n);
	auto B = Matrix(p, n);

	// you can probably turn this into a test with a pair of 1x1 matrices
	A(0,0) = +4.369503870e-07;
	A(1,0) = +1.496136406e-06;
	B(0,0) = -7.727422714e-01;
	B(1,0) = +6.347199082e-01;

	auto caller = ggqrcs::Caller<Number>(m, n, p);

	caller.A = A;
	caller.B = B;

	auto ret = caller();

	check_results(ret, A, B, caller);

	BOOST_REQUIRE_EQUAL( caller.rank, 1 );
	BOOST_REQUIRE( !caller.swapped_p );

	// computed with 2xSVD in double precision
	auto cos = Real{1.5586372e-06};

	// should fail
	BOOST_CHECK_LE(std::abs(cos - caller.alpha(0)), eps*cos);
	BOOST_CHECK_LE(std::abs(Real{1} - caller.beta(0)), eps);

	auto X = copy_X(caller);
	auto ds = ggqrcs::assemble_diagonals_like(
		Number{}, m, p, r, caller.swapped_p, caller.alpha, caller.beta);
	auto& D1 = ds.first;
	auto& D2 = ds.second;
	auto almost_A = ggqrcs::assemble_matrix(caller.U1, D1, X);
	auto almost_B = ggqrcs::assemble_matrix(caller.U2, D2, X);
	auto tol = [] (const Matrix& a) {
		return std::max(a.size1(), a.size2()) * ublas::norm_frobenius(a) * eps;
	};

	// should fail
	BOOST_CHECK_LE(ublas::norm_frobenius(A - almost_A), tol(A));
	BOOST_CHECK_LE(ublas::norm_frobenius(B - almost_B), tol(B));
}



template<
	typename Number,
	typename std::enable_if<
		std::is_fundamental<Number>::value, int
	>::type* = nullptr
>
void xGGQRCS_test_zero_dimensions_impl(
	std::size_t m, std::size_t n, std::size_t p) {
	using Storage = ublas::column_major;
	using Matrix = ublas::matrix<Number, Storage>;
	using Real = typename tools::real_from<Number>::type;
	using Integer = lapack::integer_t;

	constexpr auto nan = tools::not_a_number<Number>::value;
	constexpr auto real_nan = tools::not_a_number<Real>::value;
	constexpr auto one = std::size_t{1};

	auto k = std::max(std::min(m + p, n), one);
	auto rank = Integer{-1};
	auto swapped_p = false;
	auto lda = std::max(m, one);
	auto A = Matrix(lda, std::max(n, one), 1);
	auto ldb = std::max(p, one);
	auto B = Matrix(ldb, std::max(n, one), 1);
	auto alpha = std::vector<Real>(k, real_nan);
	auto beta = std::vector<Real>(k, real_nan);
	auto ldu1 = std::max(m, one);
	auto U1 = Matrix(ldu1, std::max(m, one), nan);
	auto ldu2 = std::max(p, one);
	auto U2 = Matrix(ldu2, std::max(p, one), nan);
	auto tol = Real{0};
	// this must be large enough not to trigger the workspace size check
	auto lwork = std::max(4 * (m + p) * n, std::size_t{128});
	auto work = std::vector<Number>(lwork, nan);
	auto iwork = std::vector<Integer>(lwork, -1);
	auto ret = lapack::xGGQRCS(
		'Y', 'Y', 'N', m, n, p, &rank, &swapped_p,
		&A(0, 0), lda, &B(0, 0), ldb,
		&alpha[0], &beta[0],
		&U1(0, 0), 1, &U2(0, 0), 1,
		&tol,
		&work[0], lwork, &iwork[0]
	);

	BOOST_REQUIRE_EQUAL(ret, 0);

	constexpr auto eps = std::numeric_limits<Real>::epsilon();
	auto nan_p = [] (const Number& x) { return tools::nan_p(x); };

	if(m > 0) {
		BOOST_REQUIRE(std::none_of(U1.data().begin(), U1.data().end(), nan_p));
		BOOST_CHECK_LE(tools::measure_isometry(U1), m * eps);
	}
	if(p > 0) {
		BOOST_REQUIRE(std::none_of(U2.data().begin(), U2.data().end(), nan_p));
		BOOST_CHECK_LE(tools::measure_isometry(U2), p * eps);
	}
}

template<
	typename Number,
	typename std::enable_if<
		!std::is_fundamental<Number>::value, int
	>::type* = nullptr
>
void xGGQRCS_test_zero_dimensions_impl(
	std::size_t m, std::size_t n, std::size_t p) {
	using Storage = ublas::column_major;
	using Matrix = ublas::matrix<Number, Storage>;
	using Real = typename tools::real_from<Number>::type;
	using Integer = lapack::integer_t;

	constexpr auto nan = tools::not_a_number<Number>::value;
	constexpr auto real_nan = tools::not_a_number<Real>::value;
	constexpr auto one = std::size_t{1};

	auto k = std::max(std::min(m + p, n), one);
	auto rank = Integer{-1};
	auto swapped_p = false;
	auto lda = std::max(m, one);
	auto A = Matrix(lda, std::max(n, one), 1);
	auto ldb = std::max(p, one);
	auto B = Matrix(ldb, std::max(n, one), 1);
	auto alpha = std::vector<Real>(k, real_nan);
	auto beta = std::vector<Real>(k, real_nan);
	auto ldu1 = std::max(m, one);
	auto U1 = Matrix(ldu1, std::max(m, one), nan);
	auto ldu2 = std::max(p, one);
	auto U2 = Matrix(ldu2, std::max(p, one), nan);
	auto tol = Real{0};
	// this must be large enough not to trigger the workspace size check
	auto lwork = std::max(4 * (m + p) * n, std::size_t{128});
	auto work = std::vector<Number>(lwork, nan);
	auto rwork = std::vector<Real>(lwork, real_nan);
	auto iwork = std::vector<Integer>(lwork, -1);
	auto ret = lapack::xGGQRCS(
		'Y', 'Y', 'N', m, n, p, &rank, &swapped_p,
		&A(0, 0), lda, &B(0, 0), ldb,
		&alpha[0], &beta[0],
		&U1(0, 0), 1, &U2(0, 0), 1,
		&tol,
		&work[0], work.size(), &rwork[0], rwork.size(), &iwork[0]
	);

	BOOST_REQUIRE_EQUAL(ret, 0);

	constexpr auto eps = std::numeric_limits<Real>::epsilon();
	auto nan_p = [] (const Number& x) { return tools::nan_p(x); };

	if(m > 0) {
		BOOST_REQUIRE(std::none_of(U1.data().begin(), U1.data().end(), nan_p));
		BOOST_CHECK_LE(tools::measure_isometry(U1), m * eps);
	}
	if(p > 0) {
		BOOST_REQUIRE(std::none_of(U2.data().begin(), U2.data().end(), nan_p));
		BOOST_CHECK_LE(tools::measure_isometry(U2), p * eps);
	}
}


BOOST_AUTO_TEST_CASE_TEMPLATE(xGGQRCS_test_zero_dimensions, Number, types)
{
	for(auto m = std::size_t{0}; m < std::size_t{2}; ++m) {
		for(auto n = std::size_t{0}; n < std::size_t{2}; ++n) {
			for(auto p = std::size_t{0}; p < std::size_t{2}; ++p) {
				xGGQRCS_test_zero_dimensions_impl<Number>(m, n, p);
			}
		}
	}
}



BOOST_AUTO_TEST_CASE_TEMPLATE(xGGQRCS_test_zero_input, Number, types)
{
	auto m = std::size_t{4};
	auto n = std::size_t{3};
	auto p = std::size_t{2};
	auto caller = ggqrcs::Caller<Number>(m, n, p);
	auto A = caller.A;
	auto B = caller.B;

	auto ret = caller();
	check_results(ret, A, B, caller);

	BOOST_CHECK_EQUAL( caller.rank, 0 );
}


BOOST_AUTO_TEST_CASE_TEMPLATE(xGGQRCS_test_rectangular_input, Number, types)
{
	for(std::size_t m : { 2, 13, 41 })
	{
		for(std::size_t n : {3, 7, 31})
		{
			for(std::size_t p : {5, 11, 17})
			{
				auto caller = ggqrcs::Caller<Number>(m, n, p);
				auto A = caller.A;
				auto B = caller.B;

				A(0,0) = 1;
				A(1,0) = 1;
				B(0,1) = 1;
				B(1,1) = 1;

				caller.A = A;
				caller.B = B;

				auto ret = caller();
				check_results(ret, A, B, caller);
			}
		}
	}
}



/**
 * This test checks the generalized singular values computed by xGGQRCS when
 * the matrices `A` and `B` differ significantly in norm.
 *
 * The GSVD allows us to decompose `A` and `B` into
 * * `A = U1 S R Q^*`,
 * * `B = U2 C R Q^*`,
 *
 * where
 * * `Q^*` is the complex-conjugate transpose of `Q`,
 * * `U1, `U2`, `Q` are unitary, and
 * * `S`, `C` are diagonal matrices which values `s_ii` and `c_ii` such that
 *   `s_ii/c_ii` is one of the generalized singular values of the matrix pencil
 *   `(A, B)`.
 *
 * To generate matrices, `A`, `B` with `A` much larger in norm than `B`, we
 * compute
 * * a random matrix `R Q^*`, and
 * * generalized singular values such that `s_ii >> c_ii`.
 */
BOOST_TEST_DECORATOR(* boost::unit_test::disabled())
BOOST_AUTO_TEST_CASE_TEMPLATE(xGGQRCS_test_singular_values, Number, types)
{
	// Numbers of the form 4n+1 are used here so that the the median as well as
	// the 25th and the 75th percentiles can be computed easily for the
	// five-number summary.
	// (Sort the values xs and use xs[0], xs[n/4], xs[n/2], xs[n/4*3], xs[n-1].)

	using Real = typename tools::real_from<Number>::type;

	auto gen = std::mt19937();

	gen.discard(1u << 17);

	std::printf(
		"%2s %8s  %44s  %44s\n",
		"d", "condition", "five-number summary norm(A)/norm(B)",
		"five-number summary relative forward error"
	);

	for(auto option = 1u; option <= 2u; ++option)
	{
		for(auto d = std::size_t{2}; d <= 65; d += 20)
		{
			auto real_nan = tools::not_a_number<Real>::value;
			auto stats_fe = ublas::vector<Real>(5, 0);
			auto cond_max = Real{0};
			auto num_iterations = std::size_t{101};
			auto rel_norms = ublas::vector<Real>(num_iterations, real_nan);

			for(auto it = std::size_t{0}; it < num_iterations; ++it)
			{
				auto m = d;
				auto n = d;
				auto p = d;
				auto r = std::min(m+p, n) - 1;
				auto k = std::min( {m, p, r, m + p - r} );

				BOOST_TEST_CONTEXT("m=" << m) {
				BOOST_TEST_CONTEXT("n=" << n) {
				BOOST_TEST_CONTEXT("p=" << p) {
				BOOST_TEST_CONTEXT("rank=" << r) {

				BOOST_VERIFY(k > 0);

				auto theta_dist = ggqrcs::ThetaDistribution<Real>(option);
				auto theta = ublas::vector<Real>(k, real_nan);

				std::generate(
					theta.begin(), theta.end(),
					[&gen, &theta_dist](){ return theta_dist(gen); }
				);
				std::sort(theta.begin(), theta.end());

				auto dummy = Number{};
				// Do not condition `X` too badly or we cannot directly compare
				// the computed generalized singular values with the generated
				// singular values.
				auto digits = std::numeric_limits<Real>::digits;
				auto cond_X = static_cast<Real>(1 << (digits/4));
				auto X = tools::make_matrix_like(dummy, r, n, cond_X, &gen);
				auto U1 = tools::make_isometric_matrix_like(dummy, m, m, &gen);
				auto U2 = tools::make_isometric_matrix_like(dummy, p, p, &gen);
				auto ds = ggqrcs::assemble_diagonals_like(dummy, m, p, r, theta);
				auto D1 = ds.first;
				auto D2 = ds.second;
				auto A = ggqrcs::assemble_matrix(U1, D1, X);
				auto B = ggqrcs::assemble_matrix(U2, D2, X);
				auto caller = ggqrcs::Caller<Number>(m, n, p);

				caller.A = A;
				caller.B = B;

				auto ret = caller();
				auto be_errors = check_results(ret, A, B, caller);
				auto norm_A = ublas::norm_frobenius(A);
				auto norm_B = ublas::norm_frobenius(B);
				auto min_be_error =
					std::min(be_errors.first/norm_A, be_errors.second/norm_B);

				BOOST_REQUIRE_LE(caller.rank, r);

				auto rank = static_cast<std::size_t>(caller.rank);
				auto l = std::min({m, p, rank, m+p-rank});
				auto alpha = ublas::subrange(caller.alpha, 0, l);
				auto beta = ublas::subrange(caller.beta, 0, l);
				auto eps = std::numeric_limits<Real>::epsilon();
				// relative forward error
				auto delta_fe = ublas::vector<Real>(l, real_nan);

				BOOST_REQUIRE_EQUAL(l, k);

				for(auto i = std::size_t{0}; i < l; ++i)
				{
					auto x = theta(i);
					auto abs = [] (Real x) { return std::abs(x); };
					auto cos = [] (Real x) { return std::cos(x); };
					auto sin = [] (Real x) { return std::sin(x); };
					auto rel_forward_error =
						(cos(x) >= sin(x)) ? abs(cos(x)-beta(i)) / cos(x) :
						                     abs(sin(x)-alpha(i)) / sin(x)
					;
					delta_fe(i) = rel_forward_error / eps;
					cond_max =
						std::max(cond_max, rel_forward_error/min_be_error);
				}

				std::sort(delta_fe.begin(), delta_fe.end());

				stats_fe(0) = std::min(stats_fe(0), delta_fe(0));
				stats_fe(1) += delta_fe(1*l/4);
				stats_fe(2) += delta_fe(2*l/4);
				stats_fe(3) += delta_fe(3*l/4);
				stats_fe(4) = std::max(stats_fe(4), delta_fe(l-1));

				rel_norms(it) = ublas::norm_frobenius(A) / ublas::norm_frobenius(B);
			}
			}
			}
			}
			}

			stats_fe(1) /= num_iterations;
			stats_fe(2) /= num_iterations;
			stats_fe(3) /= num_iterations;

			std::sort(rel_norms.begin(), rel_norms.end());

			auto rel_norm_0 =   rel_norms(num_iterations/4*0);
			auto rel_norm_25 =  rel_norms(num_iterations/4*1);
			auto rel_norm_50 =  rel_norms(num_iterations/4*2);
			auto rel_norm_75 =  rel_norms(num_iterations/4*3);
			auto rel_norm_100 = rel_norms(num_iterations/4*4);

			std::printf(
				"%2zu  %8.2e  %8.2e %8.2e %8.2e %8.2e %8.2e  %8.2e %8.2e %8.2e %8.2e %8.2e\n",
				d,
				cond_max,
				rel_norm_0, rel_norm_25, rel_norm_50, rel_norm_75, rel_norm_100,
				stats_fe(0), stats_fe(1), stats_fe(2), stats_fe(3), stats_fe(4)
			);
		}
	}
}



BOOST_TEST_DECORATOR(* boost::unit_test::disabled())
BOOST_AUTO_TEST_CASE_TEMPLATE(xGGQRCS_test_row_scaling, Number, types)
{
	using Real = typename tools::real_from<Number>::type;

	auto gen = std::mt19937();

	gen.discard(1u << 17);

	std::printf(
		"%2s %44s %44s\n",
		"d",
		"five-number summary relative backward error A",
		"five-number summary relative backward error B"
	);

	for(auto d = std::size_t{5}; d <= 45; d += 10)
	{
		auto num_iterations = std::size_t{101};
		auto real_nan = tools::not_a_number<Real>::value;
		// five-point summary backward error A
		auto stats_be_a = ublas::vector<Real>(num_iterations, real_nan);
		// five-point summary backward error B
		auto stats_be_b = ublas::vector<Real>(num_iterations, real_nan);

		for(auto it = std::size_t{0}; it < num_iterations; ++it)
		{
			auto m = d;
			auto n = d;
			auto p = 2*d;

			BOOST_TEST_CONTEXT("m=" << m) {
			BOOST_TEST_CONTEXT("n=" << n) {
			BOOST_TEST_CONTEXT("p=" << p) {
			BOOST_TEST_CONTEXT("iteration=" << it) {

			auto dummy = Number{};
			auto digits = std::numeric_limits<Real>::digits;
			auto cond = static_cast<Real>(1 << (digits/4));
			auto A = tools::make_matrix_like(dummy, m, n, cond, &gen);
			auto B = tools::make_matrix_like(dummy, p, n, cond, &gen);

			// row scaling
			auto row_norms_a = ublas::vector<Real>(m, Real{0});
			auto row_norms_b = ublas::vector<Real>(p, Real{0});

			for(auto j = std::size_t{0}; j < n; ++j)
			{
				for(auto i = std::size_t{0}; i < m; ++i)
					row_norms_a(i) = std::max(row_norms_a(i), std::abs(A(i,j)));
				for(auto i = std::size_t{0}; i < p; ++i)
					row_norms_b(i) = std::max(row_norms_b(i), std::abs(B(i,j)));
			}

			auto kappa = std::ldexp(Real{1}, digits/2);
			for(auto j = std::size_t{0}; j < n; ++j)
			{
				//for(auto i = std::size_t{0}; i < m; ++i)
				//	A(i,j) *= kappa * i / (m-1) / row_norms_a(i);
				for(auto i = std::size_t{0}; i < p; ++i)
					B(i,j) *= kappa * i / (p-1) / row_norms_b(i);
			}

			auto caller = ggqrcs::Caller<Number>(m, n, p);

			caller.A = A;
			caller.B = B;

			auto ret = caller();
			auto be_errors = check_results(ret, A, B, caller);

			stats_be_a(it) = be_errors.first;
			stats_be_b(it) = be_errors.second;
		}
		}
		}
		}
		}

		auto k = num_iterations - 1;

		std::sort(stats_be_a.begin(), stats_be_a.end());
		std::sort(stats_be_b.begin(), stats_be_b.end());
		std::printf(
			"%2zu  %8.2e %8.2e %8.2e %8.2e %8.2e  %8.2e %8.2e %8.2e %8.2e %8.2e\n",
			d,
			stats_be_a(k*0), stats_be_a(k/4), stats_be_a(k/2), stats_be_a(k/4*3), stats_be_a(k),
			stats_be_b(k*0), stats_be_b(k/4), stats_be_b(k/2), stats_be_b(k/4*3), stats_be_b(k)
		);
	}
}





template<typename Number>
void xGGQRCS_test_random_impl(
	Number dummy,
	std::size_t m, std::size_t n, std::size_t p, std::size_t r,
	std::uint64_t seed)
{
	using Real = typename tools::real_from<Number>::type;

	constexpr auto real_nan = tools::not_a_number<Real>::value;

	BOOST_TEST_CONTEXT("m=" << m) {
	BOOST_TEST_CONTEXT("n=" << n) {
	BOOST_TEST_CONTEXT("p=" << p) {
	BOOST_TEST_CONTEXT("rank=" << r) {
	BOOST_TEST_CONTEXT("seed=" << seed) {

	auto gen = std::mt19937(seed);
	auto option_dist = std::uniform_int_distribution<unsigned>(0, 2);

	gen.discard(1u << 17);

	auto option = option_dist(gen);
	auto theta_dist = ggqrcs::ThetaDistribution<Real>(option);
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

	// initialize caller
	auto ldx = m + 11;
	auto ldy = p + 5;
	auto ldu1 = m + 13;
	auto ldu2 = p + 7;
	auto caller = ggqrcs::Caller<Number>(m, n, p, ldx, ldy, ldu1, ldu2);

	ublas::subrange(caller.A, 0, m, 0, n) = A;
	ublas::subrange(caller.B, 0, p, 0, n) = B;

	auto ret = caller();

	check_results(ret, A, B, caller);

	BOOST_CHECK_LE( caller.rank, r );
}
}
}
}
}
}


BOOST_AUTO_TEST_CASE_TEMPLATE(xGGQRCS_test_random, Number, types)
{
	constexpr std::size_t dimensions[] = { 1, 2, 3, 4, 10, 20 };

	auto gen = std::mt19937();
	auto seed_dist = std::uniform_int_distribution<std::uint64_t>(0);

	gen.discard(1u << 17);

	for(auto m : dimensions)
	{
		for(auto n : dimensions)
		{
			for(auto p : dimensions)
			{
				auto max_rank = std::min(m+p, n);
				for(auto rank = std::size_t{0}; rank <= max_rank; ++rank)
				{
					for(auto iteration = 0u; iteration < 10u; ++iteration)
					{
						auto seed = seed_dist(gen);

						xGGQRCS_test_random_impl(Number{0}, m, n, p, rank, seed);
					}
				}
			}
		}
	}
}

BOOST_TEST_DECORATOR(* boost::unit_test::disabled())
BOOST_AUTO_TEST_CASE_TEMPLATE(
	xGGQRCS_test_random_infinite, Number, types)
{
	constexpr auto min_dimension = std::size_t{1};
	constexpr auto max_dimension = std::size_t{1000};

	auto master_seed = std::uintmax_t(std::time(nullptr));

	std::printf("xGGQRCS_test_random_infinite master_seed=%ju\n", master_seed);

	auto gen = std::mt19937(master_seed);
	auto dim_dist =
		std::uniform_int_distribution<std::size_t>(min_dimension,max_dimension);
	auto seed_dist = std::uniform_int_distribution<std::uint64_t>(0);

	gen.discard(1u << 17);

	auto start_time_sec = std::time(nullptr);
	auto last_time_sec = std::time_t{0};
	auto iteration = std::uintmax_t{0};

	constexpr char FMT[] = "%7jd %13ju  %3zu %3zu %3zu %4zu  %20ju\n";
	std::printf(
		"%7s %13s  %3s %3s %3s %4s  %20s\n",
		"time(s)", "iteration", "m", "n", "p", "rank", "seed"
	);

	while(true)
	{
		auto m = (dim_dist(gen) + 1) / 2;
		auto n = dim_dist(gen);
		auto p = (dim_dist(gen) + 1) / 2;
		auto max_rank = (m + p <= n) ? m + p : n;

		for(auto rank = std::size_t{0}; rank <= max_rank; ++rank, ++iteration)
		{
			auto seed = seed_dist(gen);
			auto now_sec = std::time(nullptr);
			auto second = std::time_t{1};

			if(last_time_sec + 60*second < now_sec)
			{
				auto time_passed_sec = std::intmax_t{now_sec - start_time_sec};

				std::printf(
					FMT, time_passed_sec, iteration, m, n, p, rank, seed
				);

				last_time_sec = now_sec;
			}

			xGGQRCS_test_random_impl(Number{0}, m, n, p, rank, seed);
		}
	}
}


BOOST_TEST_DECORATOR(* boost::unit_test::disabled())
BOOST_AUTO_TEST_CASE_TEMPLATE(xGGQRCS_test_xGGSVD3_comparison, Number, types)
{
	using Real = typename tools::real_from<Number>::type;
	using Matrix = ublas::matrix<Number, ublas::column_major>;

	auto master_seed = std::uintmax_t(std::time(nullptr));

	std::printf(
		"xGGQRCS_test_xGGSVD3_comparison master_seed=%ju\n",
		master_seed
	);

	auto gen = std::mt19937(master_seed);
	auto option_dist = std::uniform_int_distribution<unsigned>(0, 1);

	gen.discard(1u << 17);

	std::printf(
		"%3s %3s %3s %4s  %17s  %17s\n",
		"m", "n", "p", "rank", "median-rel-error", "max-rel-error"
	);

	for(auto dim = std::size_t{10}; dim <= 50; dim += 10)
	{
		// make odd for easy median computation
		auto num_iterations = std::size_t{1001};
		auto m = dim;
		auto n = dim;
		auto p = dim;
		auto r = std::min( m+p, n );

		BOOST_TEST_CONTEXT("m=" << m) {
		BOOST_TEST_CONTEXT("n=" << n) {
		BOOST_TEST_CONTEXT("p=" << p) {
		BOOST_TEST_CONTEXT("rank=" << r) {

		auto nan = tools::not_a_number<Number>::value;
		auto real_nan = tools::not_a_number<Real>::value;
		auto eps = std::numeric_limits<Real>::epsilon();
		auto dummy = nan;
		auto A = Matrix(1, 1, nan);
		auto B = Matrix(1, 1, nan);
		auto norm_A = real_nan;
		auto norm_B = real_nan;
		auto delta_A_qrcs = ublas::vector<Real>(num_iterations, real_nan);
		auto delta_B_qrcs = ublas::vector<Real>(num_iterations, real_nan);
		auto delta_cs_qrcs = ublas::vector<Real>(num_iterations, Real{0});
		auto delta_A_svd3 = ublas::vector<Real>(num_iterations, real_nan);
		auto delta_B_svd3 = ublas::vector<Real>(num_iterations, real_nan);
		auto delta_cs_svd3 = ublas::vector<Real>(num_iterations, Real{0});

		for(auto it = std::size_t{0}; it < num_iterations; ++it)
		{
			// set up matrices
			{
				auto k = std::min( {m, p, r, m + p - r} );
				auto theta_dist =
					std::uniform_real_distribution<Real>(0, M_PI/2);
				auto theta = ublas::vector<Real>(k, real_nan);

				std::generate(
					theta.begin(), theta.end(),
					[&gen, &theta_dist](){ return theta_dist(gen); }
				);

				auto max_log_cond_R =
					static_cast<Real>(std::numeric_limits<Real>::digits/4);
				auto cond_R = std::pow(Real{2}, max_log_cond_R);
				auto R_Qt = tools::make_matrix_like(dummy, r, n, cond_R, &gen);
				auto U1 = tools::make_isometric_matrix_like(dummy, m, m, &gen);
				auto U2 = tools::make_isometric_matrix_like(dummy, p, p, &gen);
				auto ds = ggqrcs::assemble_diagonals_like(dummy, m, p, r, theta);
				auto D1 = ds.first;
				auto D2 = ds.second;
				auto option = option_dist(gen);
				auto d = std::numeric_limits<Real>::digits - 1;
				auto w =
					(option == 0) ? std::ldexp(Real{1}, +d/2) :
					(option == 1) ? std::ldexp(Real{1}, -d/2) : real_nan
				;

				A = ggqrcs::assemble_matrix(U1, D1, R_Qt);
				B = w * ggqrcs::assemble_matrix(U2, D2, R_Qt);

				norm_A = ublas::norm_frobenius(A);
				norm_B = ublas::norm_frobenius(B);
			}

			{
				auto qrcs = ggqrcs::Caller<Number>(m, n, p);

				qrcs.A = A; qrcs.B = B;

				auto ret = qrcs();

				BOOST_VERIFY(ret == 0);

				auto X = copy_X(qrcs);
				auto ds = ggqrcs::assemble_diagonals_like(
					dummy, m, p, qrcs.rank, qrcs.swapped_p, qrcs.alpha, qrcs.beta
				);
				auto& D1 = ds.first;
				auto& D2 = ds.second;
				auto almost_A = ggqrcs::assemble_matrix(qrcs.U1, D1, X);
				auto almost_B = ggqrcs::assemble_matrix(qrcs.U2, D2, X);

				delta_A_qrcs[it] =
					ublas::norm_frobenius(A-almost_A) / (eps * norm_A);
				delta_B_qrcs[it] =
					ublas::norm_frobenius(B-almost_B) / (eps * norm_B);
			}

			{
				auto svd3 = ggsvd3::Caller<Number>(p, n, m);

				svd3.A = B; svd3.B = A;

				auto ret = svd3();

				BOOST_VERIFY(  ret == 0 );

				auto R = ggsvd3::assemble_R(svd3.k, svd3.l, svd3.A, svd3.B);
				auto ds = ggsvd3::assemble_diagonals_like(
					Number{}, p, m, svd3.k, svd3.l, svd3.alpha, svd3.beta
				);
				auto& D1 = ds.first;
				auto& D2 = ds.second;
				auto Qt = Matrix(ublas::herm(svd3.Q));
				auto almost_A = ggqrcs::assemble_matrix(svd3.U2, D2, R, Qt);
				auto almost_B = ggqrcs::assemble_matrix(svd3.U1, D1, R, Qt);

				delta_A_svd3[it] =
					ublas::norm_frobenius(A-almost_A) / (eps * norm_A);
				delta_B_svd3[it] =
					ublas::norm_frobenius(B-almost_B) / (eps * norm_B);
			}
		}

		std::sort(delta_A_qrcs.begin(), delta_A_qrcs.end());
		std::sort(delta_B_qrcs.begin(), delta_B_qrcs.end());
		std::sort(delta_A_svd3.begin(), delta_A_svd3.end());
		std::sort(delta_B_svd3.begin(), delta_B_svd3.end());

		auto k = num_iterations - 1;

		std::printf(
			"%3zu %3zu %3zu %4zu  %8.2e %8.2e  %8.2e %8.2e\n",
			m, n, p, r,
			delta_A_qrcs[k/2] / delta_A_svd3[k/2],
			delta_B_qrcs[k/2] / delta_B_svd3[k/2],
			delta_A_qrcs[k] / delta_A_svd3[k],
			delta_B_qrcs[k] / delta_B_svd3[k]
		);
	}
	}
	}
	}
	}
}


