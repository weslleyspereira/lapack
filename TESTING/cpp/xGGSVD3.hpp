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
#ifndef LAPACK_TESTS_xGGSVD3_HPP
#define LAPACK_TESTS_xGGSVD3_HPP

#include "lapack.hpp"
#include "tools.hpp"

#include <algorithm>
#include <boost/assert.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <cassert>
#include <cstdint>
#include <limits>
#include <random>
#include <type_traits>
#include <utility>


namespace lapack {
namespace ggsvd3 {

namespace ublas = boost::numeric::ublas;

using Integer = integer_t;

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
	bool compute_q_p = true;
	Integer k = -1;
	Integer l = -1;
	std::size_t m, n, p;
	std::size_t lda, ldb, ldu1, ldu2, ldqt;
	Matrix A, B;
	Matrix U1, U2, Q;
	Vector<Real> alpha;
	Vector<Real> beta;
	Vector<Number> work;
	Vector<Integer> iwork;


	Caller(std::size_t m_, std::size_t n_, std::size_t p_)
		: Caller(m_, n_, p_, m_, p_, m_, p_, n_, true, true, true)
	{}

	Caller(
		std::size_t m_, std::size_t n_, std::size_t p_,
		bool compute_u1_p_, bool compute_u2_p_, bool compute_q_p_
	) : Caller(
			m_, n_, p_, m_, p_, m_, p_, n_,
			compute_u1_p_, compute_u2_p_, compute_q_p_
		)
	{}

	Caller(
		std::size_t m_, std::size_t n_, std::size_t p_,
		std::size_t ldx_, std::size_t ldy_,
		std::size_t ldu1_, std::size_t ldu2_, std::size_t ldqt_,
		bool compute_u1_p_, bool compute_u2_p_, bool compute_q_p_
	) :
		compute_u1_p(compute_u1_p_),
		compute_u2_p(compute_u2_p_),
		compute_q_p(compute_q_p_),
		m(m_),
		n(n_),
		p(p_),
		lda(ldx_), ldb(ldy_),
		ldu1(ldu1_), ldu2(ldu2_), ldqt(ldqt_),
		A(lda, n, 0),
		B(ldb, n, 0),
		U1(ldu1, m, tools::not_a_number<Number>::value),
		U2(ldu2, p, tools::not_a_number<Number>::value),
		Q(ldqt, n, tools::not_a_number<Number>::value),
		alpha(n, tools::not_a_number<Real>::value),
		beta(n, tools::not_a_number<Real>::value),
		iwork(n, -1)
	{
		BOOST_VERIFY( m > 0 );
		BOOST_VERIFY( n > 0 );
		BOOST_VERIFY( p > 0 );
		BOOST_VERIFY( lda >= m );
		BOOST_VERIFY( ldb >= p );
		BOOST_VERIFY( ldu1 >= m );
		BOOST_VERIFY( ldu2 >= p );
		BOOST_VERIFY( ldqt >= n );

		auto nan = tools::not_a_number<Number>::value;

		// query workspace size
		auto jobu1 = compute_u1_p ? 'U' : 'N';
		auto jobu2 = compute_u2_p ? 'V' : 'N';
		auto jobq = compute_q_p ? 'Q' : 'N';
		auto lwork_opt_f = nan;
		auto ret = lapack::xGGSVD3(
			jobu1, jobu2, jobq, m, n, p, &k, &l,
			&A(0, 0), lda, &B(0, 0), ldb,
			&alpha(0), &beta(0),
			&U1(0, 0), ldu1, &U2(0, 0), ldu2, &Q(0, 0), ldqt,
			&lwork_opt_f, -1,
			&iwork(0)
		);
		BOOST_VERIFY(ret == 0);

		// resize workspace accordingly
		auto lwork_opt = static_cast<std::size_t>(std::real(lwork_opt_f));

		work.resize( lwork_opt );
		std::fill( work.begin(), work.end(), nan );
	}


	Integer operator() ()
	{
		auto jobu1 = compute_u1_p ? 'U' : 'N';
		auto jobu2 = compute_u2_p ? 'V' : 'N';
		auto jobq = compute_q_p ? 'Q' : 'N';
		return lapack::xGGSVD3(
			jobu1, jobu2, jobq, m, n, p, &k, &l,
			&A(0, 0), lda, &B(0, 0), ldb,
			&alpha(0), &beta(0),
			&U1(0, 0), ldu1, &U2(0, 0), ldu2, &Q(0, 0), ldqt,
			&work(0), work.size(),
			&iwork(0)
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
	bool compute_q_p = true;
	Integer k = -1;
	Integer l = -1;
	std::size_t m, n, p;
	std::size_t lda, ldb, ldu1, ldu2, ldqt;
	Matrix A, B;
	Matrix U1, U2, Q;
	Vector<Real> alpha, beta;
	Vector<Number> work;
	Vector<Real> rwork;
	Vector<Integer> iwork;


	Caller(std::size_t m_, std::size_t n_, std::size_t p_)
		: Caller(m_, n_, p_, m_, p_, m_, p_, n_, true, true, true)
	{}

	Caller(
		std::size_t m_, std::size_t n_, std::size_t p_,
		bool compute_u1_p_, bool compute_u2_p_, bool compute_q_p_
	) : Caller(
			m_, n_, p_, m_, p_, m_, p_, n_,
			compute_u1_p_, compute_u2_p_, compute_q_p_
		)
	{}

	Caller(
		std::size_t m_, std::size_t n_, std::size_t p_,
		std::size_t ldx_, std::size_t ldy_,
		std::size_t ldu1_, std::size_t ldu2_, std::size_t ldqt_,
		bool compute_u1_p_, bool compute_u2_p_, bool compute_q_p_
	) :
		compute_u1_p(compute_u1_p_),
		compute_u2_p(compute_u2_p_),
		compute_q_p(compute_q_p_),
		m(m_),
		n(n_),
		p(p_),
		lda(ldx_), ldb(ldy_),
		ldu1(ldu1_), ldu2(ldu2_), ldqt(ldqt_),
		A(lda, n, 0),
		B(ldb, n, 0),
		U1(ldu1, m, tools::not_a_number<Number>::value),
		U2(ldu2, p, tools::not_a_number<Number>::value),
		Q(ldqt, n, tools::not_a_number<Number>::value),
		alpha(n, tools::not_a_number<Real>::value),
		beta(n, tools::not_a_number<Real>::value),
		rwork(2*n, tools::not_a_number<Real>::value),
		iwork(n, -1)
	{
		BOOST_VERIFY( m > 0 );
		BOOST_VERIFY( n > 0 );
		BOOST_VERIFY( p > 0 );
		BOOST_VERIFY( lda >= m );
		BOOST_VERIFY( ldb >= p );
		BOOST_VERIFY( ldu1 >= m );
		BOOST_VERIFY( ldu2 >= p );
		BOOST_VERIFY( ldqt >= n );

		// query workspace size
		auto jobu1 = compute_u1_p ? 'U' : 'N';
		auto jobu2 = compute_u2_p ? 'V' : 'N';
		auto jobq = compute_q_p ? 'Q' : 'N';
		auto nan = tools::not_a_number<Number>::value;
		auto lwork_opt_f = nan;
		auto ret = lapack::xGGSVD3(
			jobu1, jobu2, jobq, m, n, p, &k, &l,
			&A(0, 0), lda, &B(0, 0), ldb,
			&alpha(0), &beta(0),
			&U1(0, 0), ldu1, &U2(0, 0), ldu2, &Q(0, 0), ldqt,
			&lwork_opt_f, -1, &rwork(0), &iwork(0) );
		BOOST_VERIFY(ret == 0);

		auto lwork_opt = static_cast<std::size_t>(std::real(lwork_opt_f));

		work.resize( lwork_opt );
		std::fill( work.begin(), work.end(), nan );
	}

	Integer operator() ()
	{
		auto jobu1 = compute_u1_p ? 'U' : 'N';
		auto jobu2 = compute_u2_p ? 'V' : 'N';
		auto jobq = compute_q_p ? 'Q' : 'N';
		return lapack::xGGSVD3(
			jobu1, jobu2, jobq, m, n, p, &k, &l,
			&A(0, 0), lda, &B(0, 0), ldb,
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
	typename Real = typename tools::real_from<Number>::type
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
	const ublas::matrix<T, Storage>& A, const ublas::matrix<T, Storage>& B)
{
	BOOST_VERIFY( A.size2() == B.size2() );

	using Matrix = ublas::matrix<T, Storage>;
	using MatrixRange = ublas::matrix_range<Matrix>;
	using ConstMatrixRange = ublas::matrix_range<const Matrix>;
	using BandedAdaptor = ublas::banded_adaptor<ConstMatrixRange>;

	auto m = A.size1();
	auto n = A.size2();
	auto r = k + l;
	auto R = Matrix(r, n, 0);

	if(r <= m)
	{
		MatrixRange R12 = ublas::subrange(R, 0, r, n-r, n);
		ConstMatrixRange X1 = ublas::subrange(A, 0, r, n-r, n);
		BandedAdaptor X1U(X1, 0, r);

		R12 = X1U;
	}
	else
	{
		MatrixRange R12 = ublas::subrange(R, 0, m, n-r, n);
		MatrixRange R22 = ublas::subrange(R, m, r, n+m-r, n);

		ConstMatrixRange X1 = ublas::subrange(A, 0, m, n-r, n);
		ConstMatrixRange Y1 = ublas::subrange(B, m-k, l, n-r+m, n);

		BandedAdaptor X1U(X1, 0, r);
		BandedAdaptor Y1U(Y1, 0, r);

		R12 = X1U;
		R22 = Y1U;
	}

	return R;
}

}
}

#endif
