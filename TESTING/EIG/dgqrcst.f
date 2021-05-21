*> \brief \b DGQRCST
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGQRCST( M, P, N, A, AF, LDA, B, BF, LDB, U, LDU, V,
*                           LDV, ALPHA, BETA, R, LDR, IWORK, WORK,
*                           LWORK, RWORK, RESULT )
*
*       .. Scalar Arguments ..
*       INTEGER            LDA, LDB, LDR, LDU, LDV, LWORK, M, N, P
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       DOUBLE PRECISION   A( LDA, * ), AF( LDA, * ), ALPHA( * ),
*      $                   B( LDB, * ), BETA( * ), BF( LDB, * ),
*      $                   R( LDR, * ), RESULT( 6 ),
*      $                   RWORK( * ), U( LDU, * ), V( LDV, * ),
*      $                   WORK( LWORK )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGQRCST tests DGGSVD3, which computes the GSVD of an M-by-N matrix A
*> and a P-by-N matrix B:
*>              A = U1 * D1 * X    and    B = U2 * D2 * X.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          P is INTEGER
*>          The number of rows of the matrix B.  P >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrices A and B.  N >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,M)
*>          The M-by-N matrix A.
*> \endverbatim
*>
*> \param[out] AF
*> \verbatim
*>          AF is DOUBLE PRECISION array, dimension (LDA,N)
*>          Details of the GSVD of A and B, as returned by DGGSVD3,
*>          see DGGSVD3 for further details.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the arrays A and AF.
*>          LDA >= max( 1,M ).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB,P)
*>          On entry, the P-by-N matrix B.
*> \endverbatim
*>
*> \param[out] BF
*> \verbatim
*>          BF is DOUBLE PRECISION array, dimension (LDB,N)
*>          Details of the GSVD of A and B, as returned by DGGSVD3,
*>          see DGGSVD3 for further details.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the arrays B and BF.
*>          LDB >= max(1,P).
*> \endverbatim
*>
*> \param[out] U
*> \verbatim
*>          U is DOUBLE PRECISION array, dimension(LDU,M)
*>          The M by M orthogonal matrix U.
*> \endverbatim
*>
*> \param[in] LDU
*> \verbatim
*>          LDU is INTEGER
*>          The leading dimension of the array U. LDU >= max(1,M).
*> \endverbatim
*>
*> \param[out] V
*> \verbatim
*>          V is DOUBLE PRECISION array, dimension(LDV,M)
*>          The P by P orthogonal matrix V.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V. LDV >= max(1,P).
*> \endverbatim
*>
*> \param[out] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION array, dimension (N)
*>
*>          The generalized singular value pairs of A and B, the
*>          ``diagonal'' matrices D1 and D2 are constructed from
*>          ALPHA and BETA, see subroutine DGGSVD3 for details.
*> \endverbatim
*>
*> \param[out] R
*> \verbatim
*>          R is DOUBLE PRECISION array, dimension(LDQ,N)
*>          The upper triangular matrix R.
*> \endverbatim
*>
*> \param[in] LDR
*> \verbatim
*>          LDR is INTEGER
*>          The leading dimension of the array R. LDR >= max(1,N).
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (N)
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (LWORK)
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK,
*>          LWORK >= max(M,P,N)*max(M,P,N).
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (max(M,P,N))
*> \endverbatim
*>
*> \param[out] RESULT
*> \verbatim
*>          RESULT is DOUBLE PRECISION array, dimension (6)
*>          The test ratios:
*>          RESULT(1) = norm( A - U1*D1*X ) / ( MAX(M,N)*norm(A)*ULP )
*>          RESULT(2) = norm( B - U2*D2*X ) / ( MAX(P,N)*norm(B)*ULP )
*>          RESULT(3) = norm( I - U'*U ) / ( M*ULP )
*>          RESULT(4) = norm( I - V'*V ) / ( P*ULP )
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup double_eig
*
*  =====================================================================
      SUBROUTINE DGQRCST( M, P, N, A, AF, LDA, B, BF, LDB, U, LDU, V,
     $                    LDV, ALPHA, BETA, R, LDR, IWORK, WORK,
     $                    LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDR, LDU, LDV, LWORK, M, N, P
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), AF( LDA, * ), ALPHA( * ),
     $                   B( LDB, * ), BETA( * ), BF( LDB, * ),
     $                   R( LDR, * ), RESULT( 4 ),
     $                   RWORK( * ), U( LDU, * ), V( LDV, * ),
     $                   WORK( LWORK )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            SWAPPED
      INTEGER            I, INFO, J, L, K, K1, K2
      DOUBLE PRECISION   ANORM, BNORM, RESID, TEMP, ULP, ULPINV, UNFL
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE, DLANSY
      EXTERNAL           DLAMCH, DLANGE, DLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DGGQRCS, DLACPY, DLASET, DSYRK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      ULP = DLAMCH( 'Precision' )
      ULPINV = ONE / ULP
      UNFL = DLAMCH( 'Safe minimum' )
*
*     Copy the matrix A to the array AF.
*
      CALL DLACPY( 'Full', M, N, A, LDA, AF, LDA )
      CALL DLACPY( 'Full', P, N, B, LDB, BF, LDB )
*
      ANORM = MAX( DLANGE( '1', M, N, A, LDA, RWORK ), UNFL )
      BNORM = MAX( DLANGE( '1', P, N, B, LDB, RWORK ), UNFL )
*
*     Factorize the matrices A and B in the arrays AF and BF.
*
      CALL DGGQRCS( 'Y', 'Y', 'Y', M, N, P, L, SWAPPED, AF, LDA,
     $              BF, LDB, ALPHA, BETA, U, LDU, V, LDV, WORK, LWORK,
     $              IWORK, INFO )
      IF ( INFO.NE.0 ) THEN
         RESULT(1) = -1
         RESULT(2) = -1
         RESULT(3) = -1
         RESULT(4) = -1
         RETURN
      ENDIF
      K  = MIN(M, P, L, M + P - L)
      K1 = MAX(L - P, 0)
      K2 = MAX(L - M, 0)
*
*     Compute A:= A - U*D1*X
*
*     X is stored in WORK(2:L*N+1)
*
      IF( .NOT.SWAPPED ) THEN
*                       k1
*     1)    A := A - [ U11 ] * [X11 X12 X13] k1
*                    [ U21 ]         
*                    [ U31 ]
         CALL DGEMM( 'No transpose', 'No transpose', M, N, K1, -ONE,
     $               U(1,1), LDU, WORK(2), L, ONE, A, LDA )
*                       k
*     2)    A := A - [ U12 ] * diag(ALPHA) * [X21 X22 X23] k
*                    [ U22 ]         
*                    [ U32 ]
         DO 20 J = 1, N
            DO 10 I = 1, K
               AF( I, J ) = ALPHA( I ) * WORK( 1 + (L*(J-1) + (K1+I)) )
   10       CONTINUE
   20    CONTINUE
*
         CALL DGEMM( 'No transpose', 'No transpose', M, N, K, -ONE,
     $               U(1,K1+1), LDU, AF, LDA, ONE, A, LDA )
      ELSE
*                       k1
*     1)    A := A - [ U13 ] * [X31 X32 X33] k1
*                    [ U23 ]         
*                    [ U33 ]
         CALL DGEMM( 'No transpose', 'No transpose', M, N, K1, -ONE,
     $               U(1,M-K1+1), LDU, WORK(2+(L-K1)), L, ONE, A, LDA )
*                       k
*     2)    A := A - [ U12 ] * diag(ALPHA) * [X21 X22 X23] k
*                    [ U22 ]         
*                    [ U32 ]
         DO 40 J = 1, N
            DO 30 I = 1, K
               AF( I, J ) = 
     $              ALPHA( I ) * WORK( 1 + (L*(J-1) + (L-K1-K+I)) )
   30       CONTINUE
   40    CONTINUE
*
         CALL DGEMM( 'No transpose', 'No transpose', M, N, K, -ONE,
     $               U(1,M-K1-K+1), LDU, AF, LDA, ONE, A, LDA )
      ENDIF
*
*     Compute norm( A - U*D1*X ) / ( MAX(1,M,N)*norm(A)*ULP ) .
*
      RESID = DLANGE( '1', M, N, A, LDA, RWORK )
*
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M, N ) ) ) / ANORM ) /
     $                 ULP
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute B:= B - V*D2*X
*
*     X is stored in WORK(2:L*N+1)
*
      IF( SWAPPED ) THEN
*                       k2
*     1)    B := B - [ V11 ] * [X11 X12 X13] k2
*                    [ V21 ]         
*                    [ V31 ]
         CALL DGEMM( 'No transpose', 'No transpose', P, N, K2, -ONE,
     $               V(1,1), LDV, WORK(2), L, ONE, B, LDB )
*                       k
*     2)    B := B - [ V12 ] * diag(BETA) * [X21 X22 X23] k
*                    [ V22 ]         
*                    [ V32 ]
         DO 60 J = 1, N
            DO 50 I = 1, K
               BF( I, J ) = BETA( I ) * WORK( 1 + (L*(J-1) + (K2+I)) )
   50       CONTINUE
   60    CONTINUE
*
         CALL DGEMM( 'No transpose', 'No transpose', P, N, K, -ONE,
     $               V(1,K2+1), LDV, BF, LDB, ONE, B, LDB )
      ELSE
*                       k2
*     1)    B := B- [ V13 ] * [X31 X32 X33] k2
*                    [ V23 ]         
*                    [ V33 ]
         CALL DGEMM( 'No transpose', 'No transpose', P, N, K2, -ONE,
     $               V(1,P-K2+1), LDV, WORK(2+(L-K2)), L, ONE, B, LDB )
*                       k
*     2)    B := B - [ V12 ] * diag(BETA) * [X21 X22 X23] k
*                    [ V22 ]         
*                    [ V32 ]
         DO 80 J = 1, N
            DO 70 I = 1, K
               BF( I, J ) = 
     $              BETA( I ) * WORK( 1 + (L*(J-1) + (L-K2-K+I)) )
   70       CONTINUE
   80    CONTINUE
*
         CALL DGEMM( 'No transpose', 'No transpose', P, N, K, -ONE,
     $               V(1,P-K2-K+1), LDV, BF, LDB, ONE, B, LDB )
      ENDIF
*
*     Compute norm( B - V*D2*X ) / ( MAX(P,N)*norm(B)*ULP ) .
*
      RESID = DLANGE( '1', P, N, B, LDB, RWORK )
      IF( BNORM.GT.ZERO ) THEN
         RESULT( 2 ) = ( ( RESID / DBLE( MAX( 1, P, N ) ) ) / BNORM ) /
     $                 ULP
      ELSE
         RESULT( 2 ) = ZERO
      END IF
*
*     Compute I - U'*U
*
      CALL DLASET( 'Full', M, M, ZERO, ONE, WORK, LDU )
      CALL DSYRK( 'Upper', 'Transpose', M, M, -ONE, U, LDU, ONE, WORK,
     $            LDU )
*
*     Compute norm( I - U'*U ) / ( M * ULP ) .
*
      RESID = DLANSY( '1', 'Upper', M, WORK, LDU, RWORK )
      RESULT( 3 ) = ( RESID / DBLE( MAX( 1, M ) ) ) / ULP
*
*     Compute I - V'*V
*
      CALL DLASET( 'Full', P, P, ZERO, ONE, WORK, LDV )
      CALL DSYRK( 'Upper', 'Transpose', P, P, -ONE, V, LDV, ONE, WORK,
     $            LDV )
*
*     Compute norm( I - V'*V ) / ( P * ULP ) .
*
      RESID = DLANSY( '1', 'Upper', P, WORK, LDV, RWORK )
      RESULT( 4 ) = ( RESID / DBLE( MAX( 1, P ) ) ) / ULP
*
      RETURN
*
*     End of DGQRCST
*
      END
