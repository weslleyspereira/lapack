/*
 * Copyright (c) 2016, 2019-2021 Christoph Conrads
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
#ifndef LAPACK_TESTS_xGGQRCS_HPP
#define LAPACK_TESTS_xGGQRCS_HPP

#include "lapack.hpp"
#include "tools.hpp"

#include <algorithm>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/test/unit_test.hpp>
#include <cassert>
#include <cstdint>
#include <limits>
#include <random>
#include <type_traits>
#include <utility>



namespace lapack {
namespace ggqrcs
{

namespace ublas = boost::numeric::ublas;

using Integer = integer_t;


/**
 * Given the matrix dimensions, the scaling factor, and the sine, cosine values
 * computed by xGGQRCS, this function assembles the GSVD diagonal matrices.
 */
template<
	typename Number,
	class Storage = ublas::column_major,
	typename Real = typename tools::real_from<Number>::type
>
std::pair< ublas::matrix<Number, Storage>, ublas::matrix<Number, Storage> >
assemble_diagonals_like(
	Number dummy,
	std::size_t m, std::size_t p, std::size_t r, bool swapped_p,
	const ublas::vector<Real>& alpha,
	const ublas::vector<Real>& beta)
{
	using Matrix = ublas::matrix<Number, Storage>;
	using IdentityMatrix = ublas::identity_matrix<Number>;

	if(swapped_p)
	{
		auto ds = assemble_diagonals_like(dummy, p, m, r, false, beta, alpha);
		return std::make_pair(ds.second, ds.first);
	}

	auto k = std::min( {m, p, r, m + p - r} );
	auto k1 = (p < r) ? r - p : std::size_t{0};
	auto k2 = (m < r) ? r - m : std::size_t{0};
	auto D1 = Matrix(m, r, 0);
	auto D2 = Matrix(p, r, 0);

	ublas::subrange(D1, 0, k1, 0, k1) = IdentityMatrix(k1);
	ublas::subrange(D2, p - k2, p, r - k2, r) = IdentityMatrix(k2);

	auto D1_22 = ublas::subrange(D1, k1, k1+k, k1, k1+k);
	auto D2_12 = ublas::subrange(D2, p-k-k2, p-k2, r-k-k2, r-k2);

	for(auto i = std::size_t{0}; i < k; ++i)
	{
		D1_22(i,i) = alpha(i);
		D2_12(i,i) = beta(i);
	}

	return std::make_pair(D1, D2);
}

/**
 * Given the matrix dimensions, the scaling factor, and the singular values in
 * radians computed by xGGQRCS, this function assembles the GSVD diagonal
 * matrices.
 *
 * The tests were written with the assumption that A = U1 S X, B = U2 C X before
 * xGGQRCS began swapping A, B as needed.
 */
template<
	typename Number,
	class Storage = ublas::column_major,
	typename Real = typename tools::real_from<Number>::type
>
std::pair< ublas::matrix<Number, Storage>, ublas::matrix<Number, Storage> >
assemble_diagonals_like(
	Number dummy,
	std::size_t m, std::size_t p, std::size_t r,
	const ublas::vector<Real>& theta)
{
	auto real_nan = tools::not_a_number<Real>::value;
	auto k = std::min( {m, p, r, m + p - r} );
	auto alpha = ublas::vector<Real>(k, real_nan);
	auto beta = ublas::vector<Real>(k, real_nan);
	auto swapped_p = true;

	for(auto i = std::size_t{0}; i < k; ++i)
	{
		alpha(i) = std::sin(theta(i));
		beta(i) = std::cos(theta(i));
	}

	return assemble_diagonals_like(dummy, m, p, r, swapped_p, alpha, beta);
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

	auto nan = tools::not_a_number<Number>::value;
	auto alpha = Number{1};
	auto beta = Number{0};
	auto DRQt = Matrix(m, n, nan);

	lapack::xGEMM(
		'N', 'N', m, n, r, alpha, &D(0,0), m, &RQt(0,0), r, beta, &DRQt(0,0), m
	);

	auto A = Matrix(m, n);

	lapack::xGEMM(
		'N', 'N', m, n, m, alpha, &U(0,0), m, &DRQt(0,0), m, beta, &A(0,0), m
	);

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

	lapack::xGEMM(
		'N', 'N', m, n, k, alpha, &R(0,0), m, &Qt(0,0), k, beta, &RQt(0,0), m
	);

	return assemble_matrix(U, D, RQt);
}



template<
	typename Number,
	class Storage,
	typename Real = typename tools::real_from<Number>::type
>
std::pair<Real, Real> check_results(
	Integer ret,
	const ublas::matrix<Number, Storage>& A,
	const ublas::matrix<Number, Storage>& B,
	Integer rank,
	bool swapped_p,
	const ublas::vector<Real> alpha,
	const ublas::vector<Real> beta,
	const ublas::matrix<Number, Storage>& U1,
	const ublas::matrix<Number, Storage>& U2,
	const ublas::matrix<Number, Storage>& X)
{
	BOOST_REQUIRE( A.size2() == B.size2() );
	BOOST_REQUIRE( A.size1() == U1.size1() );
	BOOST_REQUIRE( B.size1() == U2.size1() );
	BOOST_REQUIRE( A.size2() == X.size2() );

	constexpr auto eps = std::numeric_limits<Real>::epsilon();
	auto m = A.size1();
	auto n = A.size2();
	auto p = B.size1();

	// check scalars
	BOOST_CHECK_EQUAL( ret, 0 );

	BOOST_CHECK_GE( rank, 0 );
	BOOST_CHECK_LE( rank, std::min(m+p, n) );

	auto r = static_cast<std::size_t>(rank);
	auto k = std::min( {m, p, r, m + p - r} );

	BOOST_REQUIRE_GE(alpha.size(), k);
	BOOST_REQUIRE_GE(beta.size(), k);


	// check for NaN
	auto nan_p = [] (const Number& x) { return tools::nan_p(x); };
	BOOST_REQUIRE( std::none_of( A.data().begin(), A.data().end(), nan_p) );
	BOOST_REQUIRE( std::none_of( B.data().begin(), B.data().end(), nan_p) );
	BOOST_REQUIRE( std::none_of( X.data().begin(), X.data().end(), nan_p) );
	BOOST_REQUIRE( std::none_of( U1.data().begin(), U1.data().end(), nan_p) );
	BOOST_REQUIRE( std::none_of( U2.data().begin(), U2.data().end(), nan_p) );
	if( k > 0 )
	{
		auto nan_p = [] (const Real& x) { return tools::nan_p(x); };
		BOOST_REQUIRE( std::none_of( &alpha(0), &alpha(0)+k, nan_p) );
		BOOST_REQUIRE( std::none_of( &beta(0), &beta(0)+k, nan_p) );
	}


	// check for infinity
	auto inf_p = [] (const Number& x) { return tools::inf_p(x); };
	BOOST_REQUIRE( std::none_of( A.data().begin(), A.data().end(), inf_p) );
	BOOST_REQUIRE( std::none_of( B.data().begin(), B.data().end(), inf_p) );
	BOOST_REQUIRE( std::none_of( X.data().begin(), X.data().end(), inf_p) );
	BOOST_REQUIRE( std::none_of( U1.data().begin(), U1.data().end(), inf_p) );
	BOOST_REQUIRE( std::none_of( U2.data().begin(), U2.data().end(), inf_p) );
	if( k > 0 )
	{
		auto inf_p = [] (const Real& x) { return tools::inf_p(x); };
		BOOST_REQUIRE( std::none_of( &alpha(0), &alpha(0)+k, inf_p) );
		BOOST_REQUIRE( std::none_of( &beta(0), &beta(0)+k, inf_p) );
	}


	// check that unitary matrices are indeed unitary
	// The bound is based on Inequality (19.13) in N. J. Higham: "Accuracy and
	// Stability of Numerical Algorithms". 2002.
	BOOST_CHECK_LE( tools::measure_isometry(U1), 4 * m*m * eps );
	BOOST_CHECK_LE( tools::measure_isometry(U2), 4 * p*p * eps );


	// check the "singular values"
	for(auto i = std::size_t{0}; i < k; ++i)
	{
		BOOST_CHECK_GE(alpha(i), -eps );
		BOOST_CHECK_LE(alpha(i), Real{1}+eps );
		BOOST_CHECK_GE(beta(i), -eps );
		BOOST_CHECK_LE(beta(i), Real{1}+eps );
	}


	// reconstruct A, B from GSVD
	using Matrix = ublas::matrix<Number, Storage>;

	auto ds = assemble_diagonals_like(Number{}, m, p, r, swapped_p,alpha, beta);
	auto& D1 = ds.first;
	auto& D2 = ds.second;
	auto almost_A = assemble_matrix(U1, D1, X);
	auto almost_B = assemble_matrix(U2, D2, X);
	auto tol = [m, n, p] (const Matrix& A) {
		return 8.25 * std::max(m + p, n) * ublas::norm_frobenius(A) * eps;
	};
	auto backward_error_A = ublas::norm_frobenius(A - almost_A);
	auto backward_error_B = ublas::norm_frobenius(B - almost_B);

	BOOST_CHECK_LE(backward_error_A, tol(A));
	BOOST_CHECK_LE(backward_error_B, tol(B));

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
struct Caller
{
	using Number = Real;
	using Matrix = ublas::matrix<Number, ublas::column_major>;
	template<typename U> using Vector = ublas::vector<U>;

	bool compute_u1_p = true;
	bool compute_u2_p = true;
	bool compute_x_p = true;
	Integer rank = -1;
	bool swapped_p = false;
	std::size_t m, n, p;
	std::size_t lda, ldb, ldu1, ldu2;
	Real tol;
	Matrix A, B;
	Matrix U1, U2;
	Vector<Number> alpha, beta;
	Vector<Number> work;
	Vector<Integer> iwork;


	Caller(std::size_t m_, std::size_t n_, std::size_t p_)
		: Caller(m_, n_, p_, m_, p_, m_, p_, true, true, true)
	{}

	Caller(
		std::size_t m_, std::size_t n_, std::size_t p_,
		std::size_t lda_, std::size_t ldb_,
		std::size_t ldu1_, std::size_t ldu2_)
		: Caller(m_, n_, p_, lda_, ldb_, ldu1_, ldu2_, true, true, true)
	{}

	Caller(
		std::size_t m_, std::size_t n_, std::size_t p_,
		bool compute_u1_p_, bool compute_u2_p_, bool compute_x_p_)
		: Caller(
			m_, n_, p_, m_, p_, m_, p_,
			compute_u1_p_, compute_u2_p_, compute_x_p_
		)
	{}

	Caller(
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
		tol(-1),
		A(lda, n, 0),
		B(ldb, n, 0),
		U1(ldu1, m, tools::not_a_number<Number>::value),
		U2(ldu2, p, tools::not_a_number<Number>::value),
		alpha(n, tools::not_a_number<Real>::value),
		beta(n, tools::not_a_number<Real>::value),
		iwork(m + n + p, -1)
	{
		BOOST_VERIFY( m > 0 );
		BOOST_VERIFY( n > 0 );
		BOOST_VERIFY( p > 0 );
		BOOST_VERIFY( lda >= m );
		BOOST_VERIFY( ldb >= p );
		BOOST_VERIFY( ldu1 >= m );
		BOOST_VERIFY( ldu2 >= p );
		BOOST_VERIFY( !std::isnan(tol) );
		BOOST_VERIFY( tol <= 1 );
	}


	Integer operator() ()
	{
		constexpr auto nan = tools::not_a_number<Number>::value;

		// query workspace size
		auto jobu1 = bool2lapackjob(compute_u1_p);
		auto jobu2 = bool2lapackjob(compute_u2_p);
		auto jobx = bool2lapackjob(compute_x_p);
		auto lwork_opt_f = nan;
		auto ret = lapack::xGGQRCS(
			jobu1, jobu2, jobx, m, n, p, &rank, &swapped_p,
			&A(0, 0), lda, &B(0, 0), ldb,
			&alpha(0), &beta(0),
			&U1(0, 0), ldu1, &U2(0, 0), ldu2,
			&tol,
			&lwork_opt_f, -1, &iwork(0) );
		BOOST_REQUIRE_EQUAL( ret, 0 );

		// resize workspace accordingly
		auto lwork_opt =
			static_cast<std::size_t>(std::real(lwork_opt_f));

		work.resize( lwork_opt );
		std::fill( work.begin(), work.end(), nan );

		return lapack::xGGQRCS(
			jobu1, jobu2, jobx, m, n, p, &rank, &swapped_p,
			&A(0, 0), lda, &B(0, 0), ldb,
			&alpha(0), &beta(0),
			&U1(0, 0), ldu1, &U2(0, 0), ldu2,
			&tol,
			&work(0), work.size(), &iwork(0)
		);
	}
};


template<typename Real>
struct Caller<std::complex<Real>>
{
	using Number = std::complex<Real>;
	using Matrix = ublas::matrix<Number, ublas::column_major>;
	template<typename U> using Vector = ublas::vector<U>;

	bool compute_u1_p = true;
	bool compute_u2_p = true;
	bool compute_x_p = true;
	Integer rank = -1;
	bool swapped_p = false;
	std::size_t m, n, p;
	std::size_t lda, ldb, ldu1, ldu2;
	Real tol;
	Matrix A, B;
	Matrix U1, U2;
	Vector<Real> alpha, beta;
	Vector<Number> work;
	Vector<Real> rwork;
	Vector<Integer> iwork;


	Caller(std::size_t m_, std::size_t n_, std::size_t p_)
		: Caller(m_, n_, p_, m_, p_, m_, p_, true, true, true)
	{}

	Caller(
		std::size_t m_, std::size_t n_, std::size_t p_,
		std::size_t lda_, std::size_t ldb_,
		std::size_t ldu1_, std::size_t ldu2_)
		: Caller(m_, n_, p_, lda_, ldb_, ldu1_, ldu2_, true, true, true)
	{}

	Caller(
		std::size_t m_, std::size_t n_, std::size_t p_,
		bool compute_u1_p_, bool compute_u2_p_, bool compute_x_p_)
		: Caller(
			m_, n_, p_, m_, p_, m_, p_,
			compute_u1_p_, compute_u2_p_, compute_x_p_
		)
	{}

	Caller(
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
		tol(-1),
		A(lda, n, 0),
		B(ldb, n, 0),
		U1(ldu1, m, tools::not_a_number<Number>::value),
		U2(ldu2, p, tools::not_a_number<Number>::value),
		alpha(n, tools::not_a_number<Real>::value),
		beta(n, tools::not_a_number<Real>::value),
		iwork(m + n + p, -1)
	{
		BOOST_VERIFY( m > 0 );
		BOOST_VERIFY( n > 0 );
		BOOST_VERIFY( p > 0 );
		BOOST_VERIFY( lda >= m );
		BOOST_VERIFY( ldb >= p );
		BOOST_VERIFY( ldu1 >= m );
		BOOST_VERIFY( ldu2 >= p );
		BOOST_VERIFY( !std::isnan(tol) );
		BOOST_VERIFY( tol <= 1 );

		auto nan = tools::not_a_number<Number>::value;
		auto real_nan = tools::not_a_number<Real>::value;

		// query workspace sizes
		auto jobu1 = bool2lapackjob(compute_u1_p);
		auto jobu2 = bool2lapackjob(compute_u2_p);
		auto jobx = bool2lapackjob(compute_x_p);
		auto lwork_opt_f = nan;
		auto lrwork_opt_f = real_nan;
		auto rank = Integer{-1};
		auto ret = lapack::xGGQRCS(
			jobu1, jobu2, jobx, m, n, p, &rank, &swapped_p,
			&A(0, 0), lda, &B(0, 0), ldb,
			&alpha(0), &beta(0),
			&U1(0, 0), ldu1, &U2(0, 0), ldu2,
			&tol,
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
		return lapack::xGGQRCS(
			jobu1, jobu2, jobx, m, n, p, &rank, &swapped_p,
			&A(0, 0), lda, &B(0, 0), ldb,
			&alpha(0), &beta(0),
			&U1(0, 0), ldu1, &U2(0, 0), ldu2,
			&tol,
			&work(0), work.size(),
			&rwork(0), rwork.size(),
			&iwork(0)
		);
	}
};


template<typename Number>
ublas::matrix<Number, ublas::column_major> copy_X(
	const Caller<Number>& caller)
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
	typename Real = typename tools::real_from<Number>::type
>
std::pair<Real, Real> check_results(
	Integer ret,
	const Matrix& A, const Matrix& B,
	const Caller<Number> caller)
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
		caller.swapped_p,
		caller.alpha,
		caller.beta,
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

}
}

#endif
