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
#ifndef LAPACK_TESTS_GGQRCS_HPP
#define LAPACK_TESTS_GGQRCS_HPP

#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include "lapack.hpp"
#include "tools.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <ctime>
#include <limits>
#include <random>
#include <type_traits>
#include <utility>



namespace ublas = boost::numeric::ublas;

using Integer = lapack::integer_t;


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


template<
	typename Number,
	class Storage = ublas::column_major,
	typename Real = typename real_from<Number>::type
>
std::pair< ublas::matrix<Number, Storage>, ublas::matrix<Number, Storage> >
assemble_diagonals_like(
	Number,
	std::size_t m, std::size_t p, std::size_t r,
	const ublas::vector<Real>& theta)
{
	using Matrix = ublas::matrix<Number, Storage>;
	using IdentityMatrix = ublas::identity_matrix<Number>;
	using MatrixRange = ublas::matrix_range<Matrix>;

	auto k = std::min( {m, p, r, m + p - r} );
	auto k1 = (p < r) ? r - p : std::size_t{0};
	auto k2 = (m < r) ? r - m : std::size_t{0};

	Matrix D1(m, r, 0);
	Matrix D2(p, r, 0);

	if(k1 > 0)
	{
		MatrixRange D1_33 = ublas::subrange(D1, m - k1, m, r - k1, r);
		D1_33 = IdentityMatrix(k1);
	}

	if(k2 > 0)
	{
		MatrixRange D2_11 = ublas::subrange(D2, 0, k2, 0, k2);
		D2_11 = IdentityMatrix(k2);
	}

	if(k > 0)
	{
		MatrixRange D1_22 = ublas::subrange(D1, m-k-k1, m-k1, r-k-k1, r-k1);
		MatrixRange D2_22 = ublas::subrange(D2, k2, k2+k, k2, k2+k);

		for(std::size_t i = 0; i < k; ++i)
		{
			D1_22(i, i) = std::sin( theta(i) );
			D2_22(i, i) = std::cos( theta(i) );
		}
	}

	return std::make_pair(D1, D2);
}



template<typename Number>
ublas::matrix<Number, ublas::column_major> assemble_matrix(
	const ublas::matrix<Number, ublas::column_major>& U,
	const ublas::matrix<Number, ublas::column_major>& D,
	const ublas::matrix<Number, ublas::column_major>& RQt)
{
	BOOST_VERIFY( U.size1() == U.size2() );
	BOOST_VERIFY( U.size2() == D.size1() );
	BOOST_VERIFY( D.size2() == RQt.size1() );

	using Matrix = ublas::matrix<Number, ublas::column_major>;

	auto m = U.size1();
	auto r = RQt.size1();
	auto n = RQt.size2();

	if(r == 0)
		return Matrix(m, n, 0);

	auto nan = not_a_number<Number>::value;
	auto alpha = Number{1};
	auto beta = Number{0};
	auto DRQt = Matrix(m, n, nan);
	auto ret = lapack::gemm(
		'N', 'N', m, n, r, alpha, &D(0,0), m, &RQt(0,0), r, beta, &DRQt(0,0), m
	);

	BOOST_VERIFY( ret == 0 );

	auto A = Matrix(m, n);
	ret = lapack::gemm(
		'N', 'N', m, n, m, alpha, &U(0,0), m, &DRQt(0,0), m, beta, &A(0,0), m
	);

	BOOST_VERIFY( ret == 0 );

	return A;
}



template<typename Number>
ublas::matrix<Number, ublas::column_major> assemble_matrix(
	const ublas::matrix<Number, ublas::column_major>& U,
	const ublas::matrix<Number, ublas::column_major>& D,
	const ublas::matrix<Number, ublas::column_major>& R,
	const ublas::matrix<Number, ublas::column_major>& Qt)
{
	using Matrix = ublas::matrix<Number, ublas::column_major>;

	BOOST_VERIFY(R.size2() == Qt.size1());

	auto m = R.size1();
	auto n = Qt.size2();
	auto k = R.size2();
	auto RQt = Matrix(m, n);
	auto alpha = Number{1};
	auto beta = Number{0};
	auto ret = lapack::gemm(
		'N', 'N', m, n, k, alpha, &R(0,0), m, &Qt(0,0), k, beta, &RQt(0,0), m
	);

	BOOST_VERIFY( ret == 0 );

	return assemble_matrix(U, D, RQt);
}



template<
	typename Number,
	class Storage,
	typename Real = typename real_from<Number>::type
>
std::pair<Real, Real> check_results(
	Integer ret,
	const ublas::matrix<Number, Storage>& A,
	const ublas::matrix<Number, Storage>& B,
	Integer rank,
	const ublas::vector<Real> theta,
	const ublas::matrix<Number, Storage>& U1,
	const ublas::matrix<Number, Storage>& U2,
	const ublas::matrix<Number, Storage>& X)
{
	BOOST_REQUIRE( A.size2() == B.size2() );
	BOOST_REQUIRE( A.size1() == U1.size1() );
	BOOST_REQUIRE( B.size1() == U2.size1() );
	BOOST_REQUIRE( A.size2() == X.size2() );

	auto m = A.size1();
	auto n = A.size2();
	auto p = B.size1();
	auto eps = std::numeric_limits<Real>::epsilon();

	// check scalars
	BOOST_CHECK_EQUAL( ret, 0 );

	BOOST_CHECK_GE( rank, 0 );
	BOOST_CHECK_LE( rank, std::min(m+p, n) );

	auto r = static_cast<std::size_t>(rank);
	auto k = std::min( {m, p, r, m + p - r} );

	BOOST_REQUIRE_GE(theta.size(), k);


	// check for NaN
	bool(*isnan)(Number) = &nan_p;
	BOOST_REQUIRE( std::none_of( A.data().begin(), A.data().end(), isnan) );
	BOOST_REQUIRE( std::none_of( B.data().begin(), B.data().end(), isnan) );
	BOOST_REQUIRE( std::none_of( X.data().begin(), X.data().end(), isnan) );
	BOOST_REQUIRE( std::none_of( U1.data().begin(), U1.data().end(), isnan) );
	BOOST_REQUIRE( std::none_of( U2.data().begin(), U2.data().end(), isnan) );
	if( k > 0 )
	{
		bool(*isnan_r)(Real) = &nan_p;
		BOOST_REQUIRE( std::none_of( &theta(0), &theta(0)+k, isnan_r) );
	}


	// check for infinity
	bool(*is_inf)(Number) = &inf_p;
	BOOST_REQUIRE( std::none_of( A.data().begin(), A.data().end(), is_inf) );
	BOOST_REQUIRE( std::none_of( B.data().begin(), B.data().end(), is_inf) );
	BOOST_REQUIRE( std::none_of( X.data().begin(), X.data().end(), is_inf) );
	BOOST_REQUIRE( std::none_of( U1.data().begin(), U1.data().end(), is_inf) );
	BOOST_REQUIRE( std::none_of( U2.data().begin(), U2.data().end(), is_inf) );
	if( k > 0 )
	{
		bool(*is_inf)(Real) = &inf_p;
		BOOST_REQUIRE( std::none_of( &theta(0), &theta(0)+k, is_inf) );
	}


	// check that unitary matrices are indeed unitary
	// The bound is based on Inequality (19.13) in N. J. Higham: "Accuracy and
	// Stability of Numerical Algorithms". 2002.
	BOOST_CHECK_LE( measure_isometry(U1), 4 * m*m * eps );
	BOOST_CHECK_LE( measure_isometry(U2), 4 * p*p * eps );


	// check the "singular values"
	BOOST_CHECK_GE( theta.size(), k );

	for(auto i = std::size_t{0}; i < k; ++i)
	{
		BOOST_CHECK_GE( theta[i], 0 );
		BOOST_CHECK_LE( theta[i], Real(M_PI/2) );

		if( i > 0 )
			BOOST_CHECK_LE( theta[i-1], theta[i] );
	}


	// reconstruct A, B from GSVD
	auto ds = assemble_diagonals_like(Number{}, m, p, r, theta);
	auto& D1 = ds.first;
	auto& D2 = ds.second;
	auto almost_A = assemble_matrix(U1, D1, X);
	auto almost_B = assemble_matrix(U2, D2, X);
	auto norm_A = ublas::norm_frobenius(A);
	auto norm_B = ublas::norm_frobenius(B);
	auto norm_AB = std::sqrt(norm_A*norm_A + norm_B*norm_B);
	// xGGQRCS is only conditionally backward stable
	auto tol = 8 * std::max(m + p, n) * norm_AB * eps;
	auto backward_error_A = ublas::norm_frobenius(A - almost_A);
	auto backward_error_B = ublas::norm_frobenius(B - almost_B);

	BOOST_CHECK_LE(backward_error_A, tol);
	BOOST_CHECK_LE(backward_error_B, tol);

	return std::make_pair(backward_error_A, backward_error_B);
}



/**
 * This function maps boolean values to LAPACK flags.
 */
char bool2lapackjob(bool p)
{
	return p ? 'Y' : 'N';
}


/**
 * This structure hides the differences between real and complex implementations
 * of xGGQRCS.
 */
template<typename Real>
struct xGGQRCS_Caller
{
	using Number = Real;
	using Matrix = ublas::matrix<Number, ublas::column_major>;
	template<typename U> using Vector = ublas::vector<U>;

	bool compute_u1_p = true;
	bool compute_u2_p = true;
	bool compute_x_p = true;
	Integer rank = -1;
	std::size_t m, n, p;
	std::size_t lda, ldb, ldu1, ldu2;
	Matrix A, B;
	Matrix U1, U2;
	Vector<Number> theta;
	Vector<Number> work;
	Vector<Integer> iwork;


	xGGQRCS_Caller(std::size_t m_, std::size_t n_, std::size_t p_)
		: xGGQRCS_Caller(m_, n_, p_, m_, p_, m_, p_, true, true, true)
	{}

	xGGQRCS_Caller(
		std::size_t m_, std::size_t n_, std::size_t p_,
		std::size_t lda_, std::size_t ldb_,
		std::size_t ldu1_, std::size_t ldu2_)
		: xGGQRCS_Caller(m_, n_, p_, lda_, ldb_, ldu1_, ldu2_, true, true, true)
	{}

	xGGQRCS_Caller(
		std::size_t m_, std::size_t n_, std::size_t p_,
		bool compute_u1_p_, bool compute_u2_p_, bool compute_x_p_)
		: xGGQRCS_Caller(
			m_, n_, p_, m_, p_, m_, p_,
			compute_u1_p_, compute_u2_p_, compute_x_p_
		)
	{}

	xGGQRCS_Caller(
		std::size_t m_, std::size_t n_, std::size_t p_,
		std::size_t lda_, std::size_t ldb_,
		std::size_t ldu1_, std::size_t ldu2_,
		bool compute_u1_p_, bool compute_u2_p_, bool compute_x_p_
	) :
		compute_u1_p(compute_u1_p_),
		compute_u2_p(compute_u2_p_),
		compute_x_p(compute_x_p_),
		m(m_),
		n(n_),
		p(p_),
		lda(lda_), ldb(ldb_),
		ldu1(ldu1_), ldu2(ldu2_),
		A(lda, n, 0),
		B(ldb, n, 0),
		U1(ldu1, m, not_a_number<Number>::value),
		U2(ldu2, p, not_a_number<Number>::value),
		theta(n, not_a_number<Real>::value),
		iwork(m + n + p, -1)
	{
		BOOST_VERIFY( m > 0 );
		BOOST_VERIFY( n > 0 );
		BOOST_VERIFY( p > 0 );
		BOOST_VERIFY( lda >= m );
		BOOST_VERIFY( ldb >= p );
		BOOST_VERIFY( ldu1 >= m );
		BOOST_VERIFY( ldu2 >= p );

		auto nan = not_a_number<Number>::value;

		// query workspace size
		auto jobu1 = bool2lapackjob(compute_u1_p);
		auto jobu2 = bool2lapackjob(compute_u2_p);
		auto jobx = bool2lapackjob(compute_x_p);
		auto lwork_opt_f = nan;
		auto rank = Integer{-1};
		auto ret = lapack::ggqrcs(
			jobu1, jobu2, jobx, m, n, p, &rank,
			&A(0, 0), lda, &B(0, 0), ldb,
			&theta(0),
			&U1(0, 0), ldu1, &U2(0, 0), ldu2,
			&lwork_opt_f, -1, &iwork(0) );
		BOOST_REQUIRE_EQUAL( ret, 0 );

		// resize workspace accordingly
		auto lwork_opt =
			static_cast<std::size_t>(std::real(lwork_opt_f));

		work.resize( lwork_opt );
		std::fill( work.begin(), work.end(), nan );
	}


	Integer operator() ()
	{
		auto jobu1 = bool2lapackjob(compute_u1_p);
		auto jobu2 = bool2lapackjob(compute_u2_p);
		auto jobx = bool2lapackjob(compute_x_p);
		return lapack::ggqrcs(
			jobu1, jobu2, jobx, m, n, p, &rank,
			&A(0, 0), lda, &B(0, 0), ldb,
			&theta(0),
			&U1(0, 0), ldu1, &U2(0, 0), ldu2,
			&work(0), work.size(), &iwork(0)
		);
	}
};


template<typename Real>
struct xGGQRCS_Caller<std::complex<Real>>
{
	using Number = std::complex<Real>;
	using Matrix = ublas::matrix<Number, ublas::column_major>;
	template<typename U> using Vector = ublas::vector<U>;

	bool compute_u1_p = true;
	bool compute_u2_p = true;
	bool compute_x_p = true;
	Integer rank = -1;
	std::size_t m, n, p;
	std::size_t lda, ldb, ldu1, ldu2;
	Matrix A, B;
	Matrix U1, U2;
	Vector<Real> theta;
	Vector<Number> work;
	Vector<Real> rwork;
	Vector<Integer> iwork;


	xGGQRCS_Caller(std::size_t m_, std::size_t n_, std::size_t p_)
		: xGGQRCS_Caller(m_, n_, p_, m_, p_, m_, p_, true, true, true)
	{}

	xGGQRCS_Caller(
		std::size_t m_, std::size_t n_, std::size_t p_,
		std::size_t lda_, std::size_t ldb_,
		std::size_t ldu1_, std::size_t ldu2_)
		: xGGQRCS_Caller(m_, n_, p_, lda_, ldb_, ldu1_, ldu2_, true, true, true)
	{}

	xGGQRCS_Caller(
		std::size_t m_, std::size_t n_, std::size_t p_,
		bool compute_u1_p_, bool compute_u2_p_, bool compute_x_p_)
		: xGGQRCS_Caller(
			m_, n_, p_, m_, p_, m_, p_,
			compute_u1_p_, compute_u2_p_, compute_x_p_
		)
	{}

	xGGQRCS_Caller(
		std::size_t m_, std::size_t n_, std::size_t p_,
		std::size_t lda_, std::size_t ldb_,
		std::size_t ldu1_, std::size_t ldu2_,
		bool compute_u1_p_, bool compute_u2_p_, bool compute_x_p_
	) :
		compute_u1_p(compute_u1_p_),
		compute_u2_p(compute_u2_p_),
		compute_x_p(compute_x_p_),
		m(m_),
		n(n_),
		p(p_),
		lda(lda_), ldb(ldb_),
		ldu1(ldu1_), ldu2(ldu2_),
		A(lda, n, 0),
		B(ldb, n, 0),
		U1(ldu1, m, not_a_number<Number>::value),
		U2(ldu2, p, not_a_number<Number>::value),
		theta(n, not_a_number<Real>::value),
		iwork(m + n + p, -1)
	{
		BOOST_VERIFY( m > 0 );
		BOOST_VERIFY( n > 0 );
		BOOST_VERIFY( p > 0 );
		BOOST_VERIFY( lda >= m );
		BOOST_VERIFY( ldb >= p );
		BOOST_VERIFY( ldu1 >= m );
		BOOST_VERIFY( ldu2 >= p );

		auto nan = not_a_number<Number>::value;
		auto real_nan = not_a_number<Real>::value;

		// query workspace sizes
		auto jobu1 = bool2lapackjob(compute_u1_p);
		auto jobu2 = bool2lapackjob(compute_u2_p);
		auto jobx = bool2lapackjob(compute_x_p);
		auto lwork_opt_f = nan;
		auto lrwork_opt_f = real_nan;
		auto rank = Integer{-1};
		auto ret = lapack::ggqrcs(
			jobu1, jobu2, jobx, m, n, p, &rank,
			&A(0, 0), lda, &B(0, 0), ldb,
			&theta(0),
			&U1(0, 0), ldu1, &U2(0, 0), ldu2,
			&lwork_opt_f, -1, &lrwork_opt_f, 1, &iwork(0) );
		BOOST_REQUIRE_EQUAL( ret, 0 );

		auto lwork_opt = static_cast<std::size_t>(std::real(lwork_opt_f));
		auto lrwork_opt = static_cast<std::size_t>(std::real(lrwork_opt_f));

		work.resize( lwork_opt );
		std::fill( work.begin(), work.end(), nan );

		rwork.resize( lrwork_opt );
		std::fill( rwork.begin(), rwork.end(), real_nan );
	}


	Integer operator() ()
	{
		auto jobu1 = bool2lapackjob(compute_u1_p);
		auto jobu2 = bool2lapackjob(compute_u2_p);
		auto jobx = bool2lapackjob(compute_x_p);
		return lapack::ggqrcs(
			jobu1, jobu2, jobx, m, n, p, &rank,
			&A(0, 0), lda, &B(0, 0), ldb,
			&theta(0),
			&U1(0, 0), ldu1, &U2(0, 0), ldu2,
			&work(0), work.size(),
			&rwork(0), rwork.size(),
			&iwork(0)
		);
	}
};


template<typename Number>
ublas::matrix<Number, ublas::column_major> copy_X(
	const xGGQRCS_Caller<Number>& caller)
{
	auto m = caller.m;
	auto n = caller.n;
	auto p = caller.p;
	auto r = static_cast<std::size_t>(caller.rank);

	BOOST_VERIFY(r <= std::min(m+p, n));
	BOOST_VERIFY(1 + r*n <= caller.work.size());

	auto X = ublas::matrix<Number, ublas::column_major>(r, n);

	std::copy(&caller.work(1), &caller.work(1+r*n), X.data().begin());

	return X;
}



template<
	typename Number,
	class Matrix,
	typename Real = typename real_from<Number>::type
>
std::pair<Real, Real> check_results(
	Integer ret,
	const Matrix& A, const Matrix& B,
	const xGGQRCS_Caller<Number> caller)
{
	auto f = [] (const Matrix& A, std::size_t m, std::size_t n) {
		BOOST_VERIFY( A.size1() >= m );
		BOOST_VERIFY( A.size2() == n );
		return Matrix(ublas::subrange(A, 0, m, 0, n));
	};

	auto m = caller.m;
	auto p = caller.p;
	auto X = copy_X(caller);
	auto U1 = f(caller.U1, m, m);
	auto U2 = f(caller.U2, p, p);

	return check_results(
		ret,
		A, B,
		caller.rank,
		caller.theta,
		U1, U2,
		X
	);
}



/**
 * This random number distribution returns the radians representation of
 * generalized singular values such that `sigma = tan(theta)`.
 */
template<typename Real>
struct ThetaDistribution
{
	static_assert(std::is_fundamental<Real>::value, "");
	static_assert(!std::is_integral<Real>::value, "");

	constexpr static Real PI_2 = Real{M_PI}/2;
	constexpr static std::size_t Q =
		std::size_t{1} << ((std::numeric_limits<Real>::digits-1)/2);
	constexpr static Real nan = std::numeric_limits<Real>::quiet_NaN();

	Real min_, max_;
	std::uniform_real_distribution<Real> dist_;


	explicit ThetaDistribution(unsigned option) :
		min_(
			(option == 0) ? Real{0} :
			(option == 1) ? PI_2 / Q * 0 :
			(option == 2) ? PI_2 / Q * (Q-1) : nan
		),
		max_(
			(option == 0) ? PI_2 :
			(option == 1) ? PI_2 / Q * 1 :
			(option == 2) ? PI_2 / Q * (Q-0) : nan
		),
		dist_(std::uniform_real_distribution<Real>(min_, max_))
	{
		BOOST_VERIFY(option < 3);
		BOOST_VERIFY(std::isfinite(min_));
		BOOST_VERIFY(std::isfinite(max_));
		assert(min_ < max_);
	}




	template<class Engine>
	Real operator() (Engine& gen)
	{
		return dist_(gen);
	}
};



BOOST_AUTO_TEST_CASE_TEMPLATE(xGGQRCS_test_simple_2x2, Number, test_types)
{
	auto m = std::size_t{2};
	auto n = std::size_t{2};
	auto p = std::size_t{2};
	auto caller = xGGQRCS_Caller<Number>(m, n, p);
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

BOOST_AUTO_TEST_CASE_TEMPLATE(xGGQRCS_test_simple_1_2_3, Number, test_types)
{
	auto m = std::size_t{1};
	auto n = std::size_t{2};
	auto p = std::size_t{3};
	auto caller = xGGQRCS_Caller<Number>(m, n, p, m, p, m, p);
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



// this test does not pass with row sorting and no matrix scaling
BOOST_AUTO_TEST_CASE(xGGQRCS_test_matrix_scaling)
{
	using Number = float;

	auto m = std::size_t{1};
	auto n = std::size_t{2};
	auto p = std::size_t{10};
	auto caller = xGGQRCS_Caller<Number>(m, n, p);
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
	using Real = typename real_from<Number>::type;
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
		auto caller = xGGQRCS_Caller<Number>(m, n, p);

		caller.A = A;
		caller.B = B;

		auto ret = caller();
		BOOST_REQUIRE_EQUAL( ret, 0 );

		auto X = copy_X(caller);
		auto ds = assemble_diagonals_like(Number{}, m, p, r, caller.theta);
		auto& D1 = ds.first;
		auto& D2 = ds.second;
		auto almost_A = assemble_matrix(caller.U1, D1, X);
		auto almost_B = assemble_matrix(caller.U2, D2, X);

		BOOST_CHECK_LE(ublas::norm_frobenius(A - almost_A), tol(A));
		// should fail
		BOOST_CHECK_LE(ublas::norm_frobenius(B - almost_B), tol(B));
	}

	// try again with norm(A) = norm(B)
	{
		auto w = std::ldexp(Real{1}, 12);
		auto caller = xGGQRCS_Caller<Number>(m, n, p);

		caller.A = A;
		caller.B = w * B;

		auto ret = caller();
		BOOST_REQUIRE_EQUAL( ret, 0 );

		auto X = copy_X(caller);
		auto ds = assemble_diagonals_like(Number{}, m, p, r, caller.theta);
		auto& D1 = ds.first;
		auto& D2 = ds.second;
		auto almost_A = assemble_matrix(caller.U1, D1, X);
		auto almost_B = assemble_matrix(caller.U2, D2, X);

		BOOST_CHECK_LE(ublas::norm_frobenius(A - almost_A), tol(A));
		BOOST_CHECK_LE(ublas::norm_frobenius(w*B - almost_B), tol(w*B));
	}
}



template<
	typename Number,
	typename std::enable_if<
		std::is_fundamental<Number>::value, int
	>::type* = nullptr
>
void xGGQRCS_test_zero_dimensions_impl(Number)
{
	using Real = typename real_from<Number>::type;

	constexpr auto nan = not_a_number<Number>::value;
	constexpr auto real_nan = not_a_number<Real>::value;

	auto rank = Integer{-1};
	auto lda = 1;
	auto A = std::vector<Number>(lda*1, nan);
	auto ldb = 1;
	auto B = std::vector<Number>(ldb*1, nan);
	auto theta = std::vector<Real>(1, real_nan);
	auto lwork = 1;
	auto work = std::vector<Number>(lwork, nan);
	auto iwork = std::vector<Integer>(1, -1);
	auto f = [&] (std::size_t m, std::size_t n, std::size_t p) {
		return lapack::ggqrcs(
			'N', 'N', 'N', m, n, p, &rank,
			&A[0], lda, &B[0], ldb,
			&theta[0],
			nullptr, 1, nullptr, 1,
			&work[0], lwork, &iwork[0]
		);
	};

	BOOST_CHECK_EQUAL( f(0, 1, 1), -4 );
	BOOST_CHECK_EQUAL( f(1, 0, 1), -5 );
	BOOST_CHECK_EQUAL( f(1, 1, 0), -6 );
}

template<
	typename Number,
	typename std::enable_if<
		!std::is_fundamental<Number>::value, int
	>::type* = nullptr
>
void xGGQRCS_test_zero_dimensions_impl(Number)
{
	using Real = typename real_from<Number>::type;

	constexpr auto nan = not_a_number<Number>::value;
	constexpr auto real_nan = not_a_number<Real>::value;

	auto rank = Integer{-1};
	auto lda = 1;
	auto A = std::vector<Number>(lda*1, nan);
	auto ldb = 1;
	auto B = std::vector<Number>(ldb*1, nan);
	auto theta = std::vector<Real>(1, real_nan);
	auto lwork = 1;
	auto work = std::vector<Number>(lwork, nan);
	auto lrwork = 1;
	auto rwork = std::vector<Real>(lrwork, real_nan);
	auto iwork = std::vector<Integer>(1, -1);
	auto f = [&] (std::size_t m, std::size_t n, std::size_t p) {
		return lapack::ggqrcs(
			'N', 'N', 'N', m, n, p, &rank,
			&A[0], lda, &B[0], ldb,
			&theta[0],
			nullptr, 1, nullptr, 1,
			&work[0], lwork, &rwork[0], lrwork, &iwork[0]
		);
	};

	BOOST_CHECK_EQUAL( f(0, 1, 1), -4 );
	BOOST_CHECK_EQUAL( f(1, 0, 1), -5 );
	BOOST_CHECK_EQUAL( f(1, 1, 0), -6 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE(xGGQRCS_test_zero_dimensions, Number, test_types)
{
	xGGQRCS_test_zero_dimensions_impl(Number{0});
}



BOOST_AUTO_TEST_CASE_TEMPLATE(xGGQRCS_test_zero_input, Number, test_types)
{
	auto m = std::size_t{4};
	auto n = std::size_t{3};
	auto p = std::size_t{2};
	auto caller = xGGQRCS_Caller<Number>(m, n, p);
	auto A = caller.A;
	auto B = caller.B;

	auto ret = caller();
	check_results(ret, A, B, caller);

	BOOST_CHECK_EQUAL( caller.rank, 0 );
}


BOOST_AUTO_TEST_CASE_TEMPLATE(xGGQRCS_test_rectangular_input, Number, test_types)
{
	for(std::size_t m : { 2, 13, 41 })
	{
		for(std::size_t n : {3, 7, 31})
		{
			for(std::size_t p : {5, 11, 17})
			{
				auto caller = xGGQRCS_Caller<Number>(m, n, p);
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
BOOST_AUTO_TEST_CASE_TEMPLATE(xGGQRCS_test_singular_values, Number, test_types)
{
	// Numbers of the form 4n+1 are used here so that the the median as well as
	// the 25th and the 75th percentiles can be computed easily for the
	// five-number summary.
	// (Sort the values xs and use xs[0], xs[n/4], xs[n/2], xs[n/4*3], xs[n-1].)

	using Real = typename real_from<Number>::type;

	auto gen = std::mt19937();

	gen.discard(1u << 17);

	for(auto option = 1u; option <= 2u; ++option)
	{
		for(auto d = std::size_t{5}; d <= 65; d += 20)
		{
			auto real_nan = not_a_number<Real>::value;
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

				auto theta_dist = ThetaDistribution<Real>(option);
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
				auto X = make_matrix_like(dummy, r, n, cond_X, &gen);
				auto U1 = make_isometric_matrix_like(dummy, m, m, &gen);
				auto U2 = make_isometric_matrix_like(dummy, p, p, &gen);
				auto ds = assemble_diagonals_like(dummy, m, p, r, theta);
				auto D1 = ds.first;
				auto D2 = ds.second;
				auto A = assemble_matrix(U1, D1, X);
				auto B = assemble_matrix(U2, D2, X);
				auto caller = xGGQRCS_Caller<Number>(m, n, p);

				caller.A = A;
				caller.B = B;

				auto ret = caller();
				auto be_errors = check_results(ret, A, B, caller);
				auto min_be_error = std::min(be_errors.first, be_errors.second);

				BOOST_REQUIRE_LE(caller.rank, r);

				auto s = static_cast<std::size_t>(caller.rank);
				auto l = std::min({m, p, s, m+p-s});
				auto iota = ublas::subrange(caller.theta, 0, l);
				auto eps = std::numeric_limits<Real>::epsilon();
				// relative forward error
				auto delta_fe = ublas::vector<Real>(l, real_nan);

				BOOST_REQUIRE_EQUAL(l, k);

				for(auto i = std::size_t{0}; i < l; ++i)
				{
					auto delta = std::abs(iota(i) - theta(i));

					delta_fe(i) = delta / (eps * theta(i));
					cond_max = std::max(cond_max, delta / min_be_error);
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
BOOST_AUTO_TEST_CASE_TEMPLATE(xGGQRCS_test_row_scaling, Number, test_types)
{
	using Real = typename real_from<Number>::type;

	auto gen = std::mt19937();

	gen.discard(1u << 17);

	for(auto d = std::size_t{5}; d <= 45; d += 10)
	{
		auto num_iterations = std::size_t{101};
		auto real_nan = not_a_number<Real>::value;
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
			auto A = make_matrix_like(dummy, m, n, cond, &gen);
			auto B = make_matrix_like(dummy, p, n, cond, &gen);

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

			auto caller = xGGQRCS_Caller<Number>(m, n, p);

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
	using Real = typename real_from<Number>::type;

	constexpr auto real_nan = not_a_number<Real>::value;

	BOOST_TEST_CONTEXT("m=" << m) {
	BOOST_TEST_CONTEXT("n=" << n) {
	BOOST_TEST_CONTEXT("p=" << p) {
	BOOST_TEST_CONTEXT("rank=" << r) {
	BOOST_TEST_CONTEXT("seed=" << seed) {

	auto gen = std::mt19937(seed);
	auto option_dist = std::uniform_int_distribution<unsigned>(0, 2);

	gen.discard(1u << 17);

	auto option = option_dist(gen);
	auto theta_dist = ThetaDistribution<Real>(option);
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
	auto X = make_matrix_like(dummy, r, n, cond_X, &gen);
	auto U1 = make_isometric_matrix_like(dummy, m, m, &gen);
	auto U2 = make_isometric_matrix_like(dummy, p, p, &gen);
	auto ds = assemble_diagonals_like(dummy, m, p, r, theta);
	auto D1 = ds.first;
	auto D2 = ds.second;
	auto A = assemble_matrix(U1, D1, X);
	auto B = assemble_matrix(U2, D2, X);

	// initialize caller
	auto ldx = m + 11;
	auto ldy = p + 5;
	auto ldu1 = m + 13;
	auto ldu2 = p + 7;
	auto caller = xGGQRCS_Caller<Number>(m, n, p, ldx, ldy, ldu1, ldu2);

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


BOOST_AUTO_TEST_CASE_TEMPLATE(xGGQRCS_test_random, Number, test_types)
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
	xGGQRCS_test_random_infinite, Number, test_types)
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

	constexpr char FMT[] = "%7jd %13ju  %3zu %3zu %3zu %4zu  %20zu\n";
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
	using Real = typename real_from<Number>::type;
	using Matrix = ublas::matrix<Number, ublas::column_major>;

	auto m = std::size_t{130}; // 160
	auto n = std::size_t{131}; // 179
	auto p = std::size_t{130}; // 150
	auto ldq = m + p;
	auto gen = std::minstd_rand();
	auto Q = make_isometric_matrix_like(Real{0}, ldq, n, &gen);

	BOOST_VERIFY( is_almost_isometric(Q) );

	auto nan = not_a_number<Number>::value;
	auto real_nan = not_a_number<Real>::value;
	auto theta = ublas::vector<Real>(n, real_nan);
	auto U1 = Matrix(m, m, nan);
	auto U2 = Matrix(p, p, nan);
	auto Qt = Matrix(n, n, nan);
	auto lwork = 32 * ldq;
	auto work = ublas::vector<Number>(lwork, nan);
	auto iwork = ublas::vector<Integer>(m+p, -1);
	auto ret = lapack::uncsd2by1(
		'Y', 'Y', 'N',
		m + p, m, n,
		&Q(0, 0), ldq, &Q(m, 0), ldq,
		&theta(0),
		&U1(0,0), m, &U2(0,0), p, &Qt(0,0), n,
		&work(0), lwork, &iwork(0)
	);

	BOOST_REQUIRE_EQUAL( ret, 0 );
	BOOST_CHECK( is_almost_isometric(U1) );
	BOOST_CHECK( is_almost_isometric(U2) );
}

template<
	typename Number,
	typename std::enable_if<
		!std::is_fundamental<Number>::value, int
	>::type* = nullptr
>
void xUNCSD2BY1_regression_20200420_impl(Number)
{
	using Real = typename real_from<Number>::type;
	using Matrix = ublas::matrix<Number, ublas::column_major>;

	auto m = std::size_t{130}; // 160
	auto n = std::size_t{131}; // 179
	auto p = std::size_t{130}; // 150
	auto ldq = m + p;
	auto gen = std::minstd_rand();
	auto Q = make_isometric_matrix_like(Number{0}, ldq, n, &gen);

	BOOST_VERIFY( is_almost_isometric(Q) );

	auto nan = not_a_number<Number>::value;
	auto real_nan = not_a_number<Real>::value;
	auto theta = ublas::vector<Real>(n, real_nan);
	auto U1 = Matrix(m, m, nan);
	auto U2 = Matrix(p, p, nan);
	auto Qt = Matrix(n, n, nan);
	auto lwork = 32 * ldq;
	auto work = ublas::vector<Number>(lwork, nan);
	auto lrwork = 32 * ldq;
	auto rwork = ublas::vector<Real>(lrwork, real_nan);
	auto iwork = ublas::vector<Integer>(m+p, -1);
	auto ret = lapack::uncsd2by1(
		'Y', 'Y', 'N',
		m + p, m, n,
		&Q(0, 0), ldq, &Q(m, 0), ldq,
		&theta(0),
		&U1(0,0), m, &U2(0,0), p, &Qt(0,0), n,
		&work(0), lwork, &rwork(0), lrwork, &iwork(0)
	);

	BOOST_REQUIRE_EQUAL( ret, 0 );
	BOOST_CHECK( is_almost_isometric(U1) );
	BOOST_CHECK( is_almost_isometric(U2) );
}

BOOST_AUTO_TEST_CASE_TEMPLATE(xUNCSD2BY1_regression_20200420, Number, test_types)
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
	using Real = typename real_from<Number>::type;

	// The workspace sizes are returned as floating-point values. Thus,
	// the dimensions should not be too large if you want to check the
	// workspace query result.
	auto m = std::size_t{(1u<<8)+1};
	auto n = std::size_t{(1u<<16)+1};
	auto p = std::size_t{(1u<<19)-1};
	auto nan = not_a_number<Number>::value;
	auto real_nan = not_a_number<Real>::value;
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
	auto ret = lapack::uncsd2by1(
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
	xUNCSD2BY1_test_workspace_query_with_lrwork, Number, test_types)
{
	xUNCSD2BY1_test_workspace_query_with_lrwork_impl(Number{});
}



/**
 * This structure hides the differences between real and complex implementations
 * of xGGQRCS.
 */
template<typename Real>
struct xGGSVD3_Caller
{
	using Number = Real;
	using Matrix = ublas::matrix<Number, ublas::column_major>;
	template<typename U> using Vector = ublas::vector<U>;

	Integer k = -1;
	Integer l = -1;
	std::size_t m, n, p;
	std::size_t ldx, ldy, ldu1, ldu2, ldqt;
	Matrix X, Y;
	Matrix U1, U2, Q;
	Vector<Real> alpha;
	Vector<Real> beta;
	Vector<Number> work;
	Vector<Integer> iwork;


	xGGSVD3_Caller(std::size_t m_, std::size_t n_, std::size_t p_)
		: xGGSVD3_Caller(m_, n_, p_, m_, p_, m_, p_, n_)
	{}


	xGGSVD3_Caller(
		std::size_t m_, std::size_t n_, std::size_t p_,
		std::size_t ldx_, std::size_t ldy_,
		std::size_t ldu1_, std::size_t ldu2_, std::size_t ldqt_
	) :
		m(m_),
		n(n_),
		p(p_),
		ldx(ldx_), ldy(ldy_),
		ldu1(ldu1_), ldu2(ldu2_), ldqt(ldqt_),
		X(ldx, n, 0),
		Y(ldy, n, 0),
		U1(ldu1, m, not_a_number<Number>::value),
		U2(ldu2, p, not_a_number<Number>::value),
		Q(ldqt, n, not_a_number<Number>::value),
		alpha(n, not_a_number<Real>::value),
		beta(n, not_a_number<Real>::value),
		iwork(n, -1)
	{
		BOOST_VERIFY( m > 0 );
		BOOST_VERIFY( n > 0 );
		BOOST_VERIFY( p > 0 );
		BOOST_VERIFY( ldx >= m );
		BOOST_VERIFY( ldy >= p );
		BOOST_VERIFY( ldu1 >= m );
		BOOST_VERIFY( ldu2 >= p );
		BOOST_VERIFY( ldqt >= n );

		auto nan = not_a_number<Number>::value;

		// query workspace size
		auto lwork_opt_f = nan;
		auto ret = lapack::ggsvd3(
			'U', 'V', 'Q', m, n, p, &k, &l,
			&X(0, 0), ldx, &Y(0, 0), ldy,
			&alpha(0), &beta(0),
			&U1(0, 0), ldu1, &U2(0, 0), ldu2, &Q(0, 0), ldqt,
			&lwork_opt_f, -1,
			&iwork(0)
		);
		BOOST_REQUIRE_EQUAL( ret, 0 );

		// resize workspace accordingly
		auto lwork_opt = static_cast<std::size_t>(std::real(lwork_opt_f));

		work.resize( lwork_opt );
		std::fill( work.begin(), work.end(), nan );
	}


	Integer operator() ()
	{
		return lapack::ggsvd3(
			'U', 'V', 'Q', m, n, p, &k, &l,
			&X(0, 0), ldx, &Y(0, 0), ldy,
			&alpha(0), &beta(0),
			&U1(0, 0), ldu1, &U2(0, 0), ldu2, &Q(0, 0), ldqt,
			&work(0), work.size(),
			&iwork(0)
		);
	}
};


template<typename Real>
struct xGGSVD3_Caller<std::complex<Real>>
{
	using Number = std::complex<Real>;
	using Matrix = ublas::matrix<Number, ublas::column_major>;
	template<typename U> using Vector = ublas::vector<U>;

	Integer k = -1;
	Integer l = -1;
	std::size_t m, n, p;
	std::size_t ldx, ldy, ldu1, ldu2, ldqt;
	Matrix X, Y;
	Matrix U1, U2, Q;
	Vector<Real> alpha, beta;
	Vector<Number> work;
	Vector<Real> rwork;
	Vector<Integer> iwork;


	xGGSVD3_Caller(std::size_t m_, std::size_t n_, std::size_t p_)
		: xGGSVD3_Caller(m_, n_, p_, m_, p_, m_, p_, n_)
	{}

	xGGSVD3_Caller(
		std::size_t m_, std::size_t n_, std::size_t p_,
		std::size_t ldx_, std::size_t ldy_,
		std::size_t ldu1_, std::size_t ldu2_, std::size_t ldqt_
	) :
		m(m_),
		n(n_),
		p(p_),
		ldx(ldx_), ldy(ldy_),
		ldu1(ldu1_), ldu2(ldu2_), ldqt(ldqt_),
		X(ldx, n, 0),
		Y(ldy, n, 0),
		U1(ldu1, m, not_a_number<Number>::value),
		U2(ldu2, p, not_a_number<Number>::value),
		Q(ldqt, n, not_a_number<Number>::value),
		alpha(n, not_a_number<Real>::value),
		beta(n, not_a_number<Real>::value),
		rwork(n, not_a_number<Real>::value),
		iwork(n, -1)
	{
		BOOST_VERIFY( m > 0 );
		BOOST_VERIFY( n > 0 );
		BOOST_VERIFY( p > 0 );
		BOOST_VERIFY( ldx >= m );
		BOOST_VERIFY( ldy >= p );
		BOOST_VERIFY( ldu1 >= m );
		BOOST_VERIFY( ldu2 >= p );
		BOOST_VERIFY( ldqt >= n );

		// query workspace size
		auto nan = not_a_number<Number>::value;
		auto lwork_opt_f = nan;
		auto ret = lapack::ggsvd3(
			'U', 'V', 'Q', m, n, p, &k, &l,
			&X(0, 0), ldx, &Y(0, 0), ldy,
			&alpha(0), &beta(0),
			&U1(0, 0), ldu1, &U2(0, 0), ldu2, &Q(0, 0), ldqt,
			&lwork_opt_f, -1, &rwork(0), &iwork(0) );
		BOOST_REQUIRE_EQUAL( ret, 0 );

		auto lwork_opt = static_cast<std::size_t>(std::real(lwork_opt_f));

		work.resize( lwork_opt );
		std::fill( work.begin(), work.end(), nan );
	}

	Integer operator() ()
	{
		return lapack::ggsvd3(
			'U', 'V', 'Q', m, n, p, &k, &l,
			&X(0, 0), ldx, &Y(0, 0), ldy,
			&alpha(0), &beta(0),
			&U1(0, 0), ldu1, &U2(0, 0), ldu2, &Q(0, 0), ldqt,
			&work(0), work.size(),
			&rwork(0),
			&iwork(0)
		);
	}
};


/**
 * This function assembles the diagonal matrices of a generalized SVD with the
 * values computed by xGGSVD3.
 */
template<
	typename Number,
	class Storage = ublas::column_major,
	typename Real = typename real_from<Number>::type
>
std::pair< ublas::matrix<Number, Storage>, ublas::matrix<Number, Storage> >
assemble_diagonals_like(
	Number,
	std::size_t m, std::size_t p,
	std::size_t k, std::size_t l,
	const ublas::vector<Real>& alpha,
	const ublas::vector<Real>& beta)
{
	using Matrix = ublas::matrix<Number, Storage>;
	using IdentityMatrix = ublas::identity_matrix<Number>;

	BOOST_VERIFY( k+l <= m+p );

	auto r = k + l;
	auto D1 = Matrix(m, r, 0);
	auto D2 = Matrix(p, r, 0);

	ublas::subrange(D1, 0, k, 0, k) = IdentityMatrix(k);

	if(m >= r)
	{
		for(auto i = std::size_t{0}; i < l; ++i)
		{
			D1(k+i,k+i) = alpha(k+i);
			D2(0+i,k+i) = beta(k+i);
		}
	}
	else
	{
		assert( k <= m );

		ublas::subrange(D2, m-k, l, m, r) = IdentityMatrix(r-m);

		for(auto i = std::size_t{0}; i < m-k; ++i)
		{
			D1(k+i,k+i) = alpha(k+i);
			D2(0+i,k+i) = beta(k+i);
		}
	}

	return std::make_pair(D1, D2);
}


template<typename T, class Storage>
ublas::matrix<T, Storage> assemble_R(
	std::size_t k, std::size_t l,
	const ublas::matrix<T, Storage>& X, const ublas::matrix<T, Storage>& Y)
{
	BOOST_VERIFY( X.size2() == Y.size2() );

	using Matrix = ublas::matrix<T, Storage>;
	using MatrixRange = ublas::matrix_range<Matrix>;
	using ConstMatrixRange = ublas::matrix_range<const Matrix>;
	using BandedAdaptor = ublas::banded_adaptor<ConstMatrixRange>;

	auto m = X.size1();
	auto n = X.size2();
	auto r = k + l;
	auto R = Matrix(r, n, 0);

	if(r <= m)
	{
		MatrixRange R12 = ublas::subrange(R, 0, r, n-r, n);
		ConstMatrixRange X1 = ublas::subrange(X, 0, r, n-r, n);
		BandedAdaptor X1U(X1, 0, r);

		R12 = X1U;
	}
	else
	{
		MatrixRange R12 = ublas::subrange(R, 0, m, n-r, n);
		MatrixRange R22 = ublas::subrange(R, m, r, n+m-r, n);

		ConstMatrixRange X1 = ublas::subrange(X, 0, m, n-r, n);
		ConstMatrixRange Y1 = ublas::subrange(Y, m-k, l, n-r+m, n);

		BandedAdaptor X1U(X1, 0, r);
		BandedAdaptor Y1U(Y1, 0, r);

		R12 = X1U;
		R22 = Y1U;
	}

	return R;
}


BOOST_TEST_DECORATOR(* boost::unit_test::disabled())
BOOST_AUTO_TEST_CASE_TEMPLATE(xGGQRCS_test_xGGSVD3_comparison, Number, test_types)
{
	using Real = typename real_from<Number>::type;
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

		auto nan = not_a_number<Number>::value;
		auto real_nan = not_a_number<Real>::value;
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
				auto R_Qt = make_matrix_like(dummy, r, n, cond_R, &gen);
				auto U1 = make_isometric_matrix_like(dummy, m, m, &gen);
				auto U2 = make_isometric_matrix_like(dummy, p, p, &gen);
				auto ds = assemble_diagonals_like(dummy, m, p, r, theta);
				auto D1 = ds.first;
				auto D2 = ds.second;
				auto option = option_dist(gen);
				auto d = std::numeric_limits<Real>::digits - 1;
				auto w =
					(option == 0) ? std::ldexp(Real{1}, +d/2) :
					(option == 1) ? std::ldexp(Real{1}, -d/2) : real_nan
				;

				A = assemble_matrix(U1, D1, R_Qt);
				B = w * assemble_matrix(U2, D2, R_Qt);

				norm_A = ublas::norm_frobenius(A);
				norm_B = ublas::norm_frobenius(B);
			}

			{
				auto qrcs = xGGQRCS_Caller<Number>(m, n, p);

				qrcs.A = A; qrcs.B = B;

				auto ret = qrcs();

				BOOST_VERIFY(ret == 0);

				auto X = copy_X(qrcs);
				auto ds = assemble_diagonals_like(
					dummy, m, p, qrcs.rank, qrcs.theta
				);
				auto& D1 = ds.first;
				auto& D2 = ds.second;
				auto almost_A = assemble_matrix(qrcs.U1, D1, X);
				auto almost_B = assemble_matrix(qrcs.U2, D2, X);

				delta_A_qrcs[it] =
					ublas::norm_frobenius(A-almost_A) / (eps * norm_A);
				delta_B_qrcs[it] =
					ublas::norm_frobenius(B-almost_B) / (eps * norm_B);
			}

			{
				auto svd3 = xGGSVD3_Caller<Number>(p, n, m);

				svd3.X = B; svd3.Y = A;

				auto ret = svd3();

				BOOST_VERIFY(  ret == 0 );

				auto R = assemble_R(svd3.k, svd3.l, svd3.X, svd3.Y);
				auto ds = assemble_diagonals_like(
					Number{}, p, m, svd3.k, svd3.l, svd3.alpha, svd3.beta
				);
				auto& D1 = ds.first;
				auto& D2 = ds.second;
				auto Qt = Matrix(ublas::herm(svd3.Q));
				auto almost_A = assemble_matrix(svd3.U2, D2, R, Qt);
				auto almost_B = assemble_matrix(svd3.U1, D1, R, Qt);

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

#endif
