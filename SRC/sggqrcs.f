*> \brief <b> SGGQRCS computes the singular value decomposition (SVD) for OTHER matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SGGQRCS + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sggqrcs.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sggqrcs.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sggqrcs.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGGQRCS( JOBU1, JOBU2, JOBX, M, N, P, RANK, SWAPPED,
*                           A, LDA, B, LDB,
*                           ALPHA, BETA,
*                           U1, LDU1, U2, LDU2,
*                           TOL,
*                           WORK, LWORK, IWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBU1, JOB2, JOBX
*       INTEGER            INFO, LDA, LDB, LDU1, LDU2, M, N, P, RANK, LWORK
*       REAL               TOL
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       REAL               A( LDA, * ), B( LDB, * ),
*      $                   ALPHA( N ), BETA( N ),
*      $                   U1( LDU1, * ), U2( LDU2, * ),
*      $                   WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGGQRCS computes the generalized singular value decomposition (GSVD)
*> of an M-by-N real matrix A and P-by-N real matrix B:
*>
*>       A = U1 * D1 * X,           B = U2 * D2 * X
*>
*> where U1 and U2 are orthogonal matrices. SGGQRCS uses the QR
*> factorization with column pivoting and the 2-by-1 CS decomposition to
*> compute the GSVD.
*>
*> Let RANK be the effective numerical rank of the matrix (A**T,B**T)**T,
*> then X is a RANK-by-N nonsingular matrix, D1 and D2 are M-by-RANK and
*> P-by-RANK "diagonal" matrices. If SWAPPED is false, then D1 and D2 are
*> of the of the following structures, respectively:
*>
*>                 K1  K
*>            K1 [ I   0   0 ]
*>       D1 = K  [ 0   C   0 ]
*>               [ 0   0   0 ]
*>
*>                     K   K2
*>               [ 0   0   0 ]
*>       D2 = K  [ 0   S   0 ]
*>            K2 [ 0   0   I ]
*>
*> where
*>
*>   K  = MIN(M, P, RANK, M + P - RANK),
*>   K1 = MAX(RANK - P, 0),
*>   K2 = MAX(RANK - M, 0),
*>   C  = diag( ALPHA(1), ..., ALPHA(K) ),
*>   S  = diag( BETA(1), ..., BETA(K) ), and
*>   C^2 + S^2 = I.
*>
*> If SWAPPED is true, then D1 and D2 are of the of the following
*> structures, respectively:
*>
*>                     K   K1
*>               [ 0   0   0 ]
*>       D1 = K  [ 0   S   0 ]
*>            K1 [ 0   0   I ]
*>
*>                 K2  K
*>            K2 [ I   0   0 ]
*>       D2 = K  [ 0   C   0 ]
*>               [ 0   0   0 ]
*>
*> where
*>
*>   S  = diag( ALPHA(1), ..., ALPHA(K) ),
*>   C  = diag( BETA(1), ..., BETA(K) ), and
*>   C^2 + S^2 = I.
*>
*> The routine computes C, S and optionally the matrices U1, U2, and X.
*> On exit, X is stored in WORK( 2:RANK*N+1 ).
*>
*> If B is an N-by-N nonsingular matrix, then the GSVD of the matrix
*> pair (A, B) implicitly gives the SVD of A*inv(B):
*>
*>       A*inv(B) = U1*(D1*inv(D2))*U2**T.
*>
*> If (A**T,B**T)**T  has orthonormal columns, then the GSVD of A and B
*> is also equal to the CS decomposition of A and B. Furthermore, the
*> GSVD can be used to derive the solution of the eigenvalue problem:
*>
*>       A**T*A x = lambda * B**T*B x.
*>
*> In some literature, the GSVD of A and B is presented in the form
*>
*>       A = U1*D1*( 0 R )*Q**T,    B = U2*D2*( 0 R )*Q**T
*>
*> where U1, U2, and Q are orthogonal matrices. This latter GSVD form is
*> computed directly by DGGSVD3. It is possible to convert between the
*> two representations by calculating the RQ decomposition of X but this
*> is not recommended for reasons of numerical stability.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOBU1
*> \verbatim
*>          JOBU1 is CHARACTER*1
*>          = 'Y':  Orthogonal matrix U1 is computed;
*>          = 'N':  U1 is not computed.
*> \endverbatim
*>
*> \param[in] JOBU2
*> \verbatim
*>          JOBU2 is CHARACTER*1
*>          = 'Y':  Orthogonal matrix U2 is computed;
*>          = 'N':  U2 is not computed.
*> \endverbatim
*>
*> \param[in] JOBX
*> \verbatim
*>          JOBX is CHARACTER*1
*>          = 'Y':  Matrix X is computed;
*>          = 'N':  X is not computed.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 1.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrices A and B.  N >= 1.
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          P is INTEGER
*>          The number of rows of the matrix B.  P >= 1.
*> \endverbatim
*>
*> \param[out] RANK
*> \verbatim
*>          RANK is INTEGER
*>          On exit, the effective numerical rank of the matrix
*>          (A**T, B**T)**T.
*> \endverbatim
*>
*> \param[out] SWAPPED
*> \verbatim
*>          RANK is LOGICAL
*>          On exit, SWAPPED is true if SGGQRCS swapped the input
*>          matrices A, B and computed the GSVD of (B, A); false
*>          otherwise.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is REAL array, dimension (LDB,N)
*>          On entry, the P-by-N matrix B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B. LDB >= max(1,P).
*> \endverbatim
*>
*> \param[out] ALPHA
*> \verbatim
*>          ALPHA is REAL array, dimension (N)
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          BETA is REAL array, dimension (N)
*>
*>          On exit, ALPHA and BETA contain the K generalized singular
*>          value pairs of A and B.
*> \endverbatim
*>
*> \param[out] U1
*> \verbatim
*>          U1 is REAL array, dimension (LDU1,M)
*>          If JOBU1 = 'Y', U1 contains the M-by-M orthogonal matrix U1.
*>          If JOBU1 = 'N', U1 is not referenced.
*> \endverbatim
*>
*> \param[in] LDU1
*> \verbatim
*>          LDU1 is INTEGER
*>          The leading dimension of the array U1. LDU1 >= max(1,M) if
*>          JOBU1 = 'Y'; LDU1 >= 1 otherwise.
*> \endverbatim
*>
*> \param[out] U2
*> \verbatim
*>          U2 is REAL array, dimension (LDU2,P)
*>          If JOBU2 = 'Y', U2 contains the P-by-P orthogonal matrix U2.
*>          If JOBU2 = 'N', U2 is not referenced.
*> \endverbatim
*>
*> \param[in] LDU2
*> \verbatim
*>          LDU2 is INTEGER
*>          The leading dimension of the array U2. LDU2 >= max(1,P) if
*>          JOBU2 = 'Y'; LDU2 >= 1 otherwise.
*> \endverbatim
*>
*> \param[in,out] TOL
*> \verbatim
*>          TOL is REAL
*>          This user-provided tolerance is used for the rank determination
*>          of the matrix G = (A**T, W*B**T)**T, see the documentation
*>          of ABSTOL for details.
*>
*>          If TOL < 0, then the tolerance will be determined
*>          automatically and this should be the default choice for most
*>          users. Otherwise, the user must provide a value in the
*>          closed interval [0, 1].
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>
*>          If LWORK = -1, then a workspace query is assumed; the
*>          routine only calculates the optimal size of the WORK array,
*>          returns this value as the first entry of the WORK array, and
*>          no error message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (M + N + P)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          > 0:  SBBCSD did not converge. For further details, see
*>                subroutine SORCSDBY1.
*> \endverbatim
*
*> \par Internal Parameters:
*  =========================
*>
*> \verbatim
*>  W       REAL
*>          W is a radix power chosen such that the Frobenius norm of A
*>          and W*B are within SQRT(RADIX) and 1/SQRT(RADIX) of each
*>          other.
*>
*>  ABSTOL  REAL
*>          Let G = (A**T, W*B**T)**T. ABSTOL is the threshold to determine
*>          the effective rank of G. Generally, it is set to
*>                   ABSTOL = TOL * MAX( M + P, N ) * norm(G),
*>          where norm(G) is the Frobenius norm of G.
*>          The size of ABSTOL may affect the size of backward error of the
*>          decomposition.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Christoph Conrads (https://christoph-conrads.name)
*
*> \ingroup realGEsing
*
*> \par Contributors:
*  ==================
*>
*>     Christoph Conrads (https://christoph-conrads.name)
*>
*
*> \par Further Details:
*  =====================
*>
*>  SGGQRCS should be significantly faster than SGGSVD3 for large
*>  matrices because the matrices A and B are reduced to a pair of
*>  well-conditioned bidiagonal matrices instead of pairs of upper
*>  triangular matrices. On the downside, SGGQRCS requires a much larger
*>  workspace whose dimension must be queried at run-time. SGGQRCS also
*>  offers no guarantees which of the two possible diagonal matrices
*>  is used for the matrix factorization.
*>
*  =====================================================================
      RECURSIVE SUBROUTINE SGGQRCS( JOBU1, JOBU2, JOBX,
     $                              HINTPREPA, HINTPREPB,
     $                              M, N, P, RANK,
     $                              SWAPPED,
     $                              A, LDA, B, LDB,
     $                              ALPHA, BETA,
     $                              U1, LDU1, U2, LDU2, X, LDX,
     $                              TOL,
     $                              WORK, LWORK,
     $                              IWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --
*
      IMPLICIT NONE
*     .. Scalar Arguments ..
      LOGICAL            SWAPPED
      CHARACTER          JOBU1, JOBU2, JOBX, HINTPREPA, HINTPREPB
      INTEGER            INFO, LDA, LDB, LDU1, LDU2, LDX,
     $                   M, N, P, RANK, LWORK
      REAL               TOL
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      REAL               A( LDA, * ), B( LDB, * ),
     $                   ALPHA( N ), BETA( N ),
     $                   U1( LDU1, * ), U2( LDU2, * ), X( LDX, * ),
     $                   WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            PREPROCESSA, PREPROCESSB,
     $                   WANTU1, WANTU2, WANTX, LQUERY
      INTEGER            I, J,
     $                   K, K1, K2, KP, K1P, K2P,
     $                   RANKMAXA, RANKMAXB, RANKMAXG, ROWSA, ROWSB,
     $                   ITAUA, ITAUB, ITAUG, IG, ISCRATCH,
     $                   IG11, IG21, IG22, LDG,
     $                   LWKMIN, LWKOPT
      REAL               BASE, ULP, UNFL,
     $                   NORMA, NORMB, NORMG,
     $                   ABSTOLA, ABSTOLB, ABSTOLG,
     $                   THETA, IOTA, W,
     $                   NAN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH, SLANGE
      EXTERNAL           LSAME, SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM, SGEQP3, SLACPY, SLAPMT, SLASCL,
     $                   SLASET, SORGQR, SORCSD2BY1, SORMQR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ACOS, COS, ISNAN, MAX, MIN, SIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     IWORK stores the column permutations computed by xGEQP3.
*     Columns J where IWORK( J ) is non-zero are permuted to the front
*     so IWORK must be set to zero before every call to xGEQP3.
*
*
*     Test the input arguments
*
*     RANKMAXG is needed in a test below
      LQUERY = LWORK.EQ.-1
      RANKMAXA = MIN( M, N )
      RANKMAXB = MIN( P, N )
      RANKMAXG = MIN( RANKMAXA + RANKMAXB, N )
*
      INFO = 0
      IF( .NOT.( LSAME( JOBU1, 'Y' ) .OR. LSAME( JOBU1, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LSAME(JOBU2, 'Y') .OR. LSAME(JOBU2, 'N') ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( LSAME(JOBX, 'Y') .OR. LSAME(JOBX, 'N') ) ) THEN
         INFO = -3
      ELSE IF( .NOT.( LSAME( HINTPREPA,'Y' )
     $                .OR. LSAME( HINTPREPA, '?' )
     $                .OR. LSAME( HINTPREPA, 'N' ) ) ) THEN
         INFO = -4
      ELSE IF( .NOT.( LSAME( HINTPREPB,'Y' )
     $                .OR. LSAME( HINTPREPB, '?' )
     $                .OR. LSAME( HINTPREPB, 'N' ) ) ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 ) THEN
         INFO = -7
      ELSE IF( P.LT.0 ) THEN
         INFO = -8
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -12
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -14
      ELSE IF( LDU1.LT.1 .OR. ( WANTU1 .AND. LDU1.LT.M ) ) THEN
         INFO = -18
      ELSE IF( LDU2.LT.1 .OR. ( WANTU2 .AND. LDU2.LT.P ) ) THEN
         INFO = -20
      ELSE IF( LDX.LT.1 .OR. ( WANTX .AND. LDX.LT.RANKMAXG ) ) THEN
         INFO = -22
      ELSE IF( ISNAN(TOL) .OR. TOL.GT.1.0E0 ) THEN
         INFO = -23
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -25
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGGQRCS', -INFO )
         RETURN
      END IF
*
*     Initialize variables
*
      SWAPPED = .FALSE.
      PREPROCESSA =
     $ M.GT.N .OR. ( M.GT.0 .AND. .NOT.LSAME(HINTPREPA, 'N') )
      PREPROCESSB =
     $ P.GT.N .OR. ( P.GT.0 .AND. .NOT.LSAME(HINTPREPB, 'N') )
      PREPROCESSA = .FALSE.
      PREPROCESSB = .FALSE.
      WANTU1 = LSAME( JOBU1, 'Y' )
      WANTU2 = LSAME( JOBU2, 'Y' )
      WANTX = LSAME( JOBX, 'Y' )
*
      RANK = -1
      IF( PREPROCESSA ) THEN
         ROWSA = MIN( M, N )
      ELSE
         ROWSA = M
      ENDIF
      IF( PREPROCESSB ) THEN
         ROWSB = MIN( P, N )
      ELSE
         ROWSB = P
      ENDIF
*     The leading dimension must never be zero
      LDG = MAX( ROWSA + ROWSB, 1 )
*     Compute offsets into workspace
      ITAUA = 1
      IF( PREPROCESSA ) THEN
         ITAUB = ITAUA + RANKMAXA
      ELSE
         ITAUB = ITAUA
      ENDIF
      IF( PREPROCESSB ) THEN
         IG = ITAUB + RANKMAXB
      ELSE
         IG = ITAUB
      ENDIF
      ITAUG = IG + LDG * N
      ISCRATCH = ITAUG + RANKMAXG
      IG11 = IG
      IG21 = IG + ROWSA
      IG22 = IG + LDG * M + M
      LWKMIN = -1
      LWKOPT = -1
*
      BASE = SLAMCH( 'B' )
      ULP = SLAMCH( 'Precision' )
      UNFL = SLAMCH( 'Safe Minimum' )
      IF( TOL.LT.0.0E0 .AND. .NOT.LQUERY ) THEN
         TOL = ULP
      ENDIF
*
      NORMA = SLANGE( 'F', M, N, A, LDA, WORK )
      NORMB = SLANGE( 'F', P, N, B, LDB, WORK )
      ABSTOLA = TOL * MAX( M, N ) * MAX( NORMA, UNFL )
      ABSTOLB = TOL * MAX( P, N ) * MAX( NORMB, UNFL )
      ABSTOLG = -1
      THETA = -1
      IOTA = -1
      W = -1
*
      IF( ISNAN(NORMA) ) THEN
         INFO = 101
         CALL XERBLA( 'SGGQRCS', INFO )
         RETURN
      ENDIF
*
      IF( ISNAN(NORMB) ) THEN
         INFO = 102
         CALL XERBLA( 'SGGQRCS', INFO )
         RETURN
      ENDIF
*
*     Make sure A is the matrix smaller in norm
*
      IF( NORMA.GT.SQRT( 2.0E0 ) * NORMB ) THEN
         CALL SGGQRCS( JOBU2, JOBU1, JOBX, HINTPREPA, HINTPREPB,
     $                 P, N, M, RANK,
     $                 SWAPPED,
     $                 B, LDB, A, LDA,
     $                 BETA, ALPHA,
     $                 U2, LDU2, U1, LDU1, X, LDX,
     $                 TOL,
     $                 WORK, LWORK, IWORK, INFO )
         SWAPPED = .TRUE.
         RETURN
      ENDIF
*
*     Past this point, we know that
*     * NORMA <= NORMB (almost)
*     * W >= 1
*     * ALPHA will contain cosine values at the end
*     * BETA will contain sine values at the end
*
*
*     Inform caller about pre-processing decisions
*
      IF( PREPROCESSA ) THEN
         HINTPREPA = 'Y'
      ELSE
         HINTPREPA = 'N'
      ENDIF
      IF( PREPROCESSB ) THEN
         HINTPREPB = 'Y'
      ELSE
         HINTPREPB = 'N'
      ENDIF
*
*     Compute workspace
*
      LWKMIN = 0
      LWKOPT = 0
*
      IF( PREPROCESSA ) THEN
         CALL SGEQP3( M, N, A, LDA, IWORK, WORK( ITAUA ),
     $                WORK, -1, INFO )
         LWKMIN = MAX( LWKMIN, 3 * N + 1 )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
      ENDIF
*
      IF( PREPROCESSB ) THEN
         CALL SGEQP3( P, N, B, LDB, IWORK, WORK( ITAUB ),
     $                WORK, -1, INFO )
         LWKMIN = MAX( LWKMIN, 3 * N + 1 )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
      ENDIF
*
      CALL SGEQP3( ROWSA + ROWSB, N, WORK( IG ), LDG, IWORK,
     $             WORK( ITAUG ), WORK, -1, INFO )
      LWKMIN = MAX( LWKMIN, 3 * N + 1 )
      LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
*
      CALL SORGQR( ROWSA + ROWSB, RANKMAXG, RANKMAXG, WORK( IG ), LDG,
     $             WORK( ITAUG ), WORK, -1, INFO )
      LWKMIN = MAX( LWKMIN, RANKMAXG )
      LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
*     Add workspace for xGGQRCS
      LWKMIN = LWKMIN + ISCRATCH
      LWKOPT = LWKOPT + ISCRATCH
*
      CALL SORCSD2BY1( JOBU1, JOBU2, JOBX,
     $                 ROWSA + ROWSB, ROWSA, RANKMAXG,
     $                 WORK( IG ), LDG, WORK( IG21 ), LDG,
     $                 BETA,
     $                 U1, LDU1, U2, LDU2, X, LDX,
     $                 WORK, -1, IWORK, INFO )
*     By the time xORCSD2BY1 is called, TAU(G) is not needed anymore
      LWKMIN = MAX( LWKMIN, INT( WORK( 1 ) ) + ITAUG )
      LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) + ITAUG )
*     Check workspace size
      IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
         INFO = -25
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGGQRCS', -INFO )
         RETURN
      END IF
      IF( LQUERY ) THEN
         WORK( 1 ) = REAL( LWKOPT )
         RETURN
      ENDIF
*
*     DEBUG
*
      NAN = 0.0E0
      NAN = 0.0E0 / NAN
      IWORK( :M+N+P ) = -1
      WORK( :LWORK ) = NAN
      ALPHA( :N ) = NAN
      BETA( :N ) = NAN
      IF( WANTU1 ) THEN
         U1( :, :M ) = NAN
      ENDIF
      IF( WANTU2 ) THEN
         U2( :, :P ) = NAN
      ENDIF
      IF( WANTX ) THEN
         X( :, 1:N ) = NAN
      ENDIF
*
*     Set scaling factor W such that norm(A) \approx norm(B)
*
      IF( NORMA.EQ.0.0E0 ) THEN
         W = 1.0E0
      ELSE
         W = BASE ** INT( LOG( NORMB / NORMA ) / LOG( BASE ) )
      END IF
*
*     Attempt to remove unnecessary matrix rows
*     Copy matrices A, B or their full-rank factors, respectively, into
*     the LDG x N matrix G
*
      IF( PREPROCESSA ) THEN
         IWORK( 1:N ) = 0
         CALL SGEQP3( M, N, A, LDA, IWORK, WORK( ITAUA ),
     $                WORK( ITAUB ), LWORK - ITAUB + 1, INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
*        Determine rank of A
         ROWSA = 0
         DO I = 1, MIN( M, N )
            IF( ABS( A( I, I ) ).LE.ABSTOLA ) THEN
               EXIT
            ENDIF
            ROWSA = ROWSA + 1
         END DO
         IG21 = IG + ROWSA
*        Scale, copy full rank part into G
         CALL SLASCL( 'U', -1, -1, 1.0E0, W, ROWSA, N, A, LDA, INFO )
         IF ( INFO.NE.0 ) THEN
            RETURN
         END IF
         CALL SLASET( 'L', ROWSA, N, 0.0E0, 0.0E0, WORK( IG11 ), LDG )
         CALL SLACPY( 'U', ROWSA, N, A, LDA, WORK( IG11 ), LDG )
         CALL SLAPMT( .FALSE., ROWSA, N, WORK( IG11 ), LDG, IWORK )
*        Initialize U1 although xORCSDB2BY1 will partially overwrite this
         IF( WANTU1 ) THEN
             CALL SLASET( 'A', M, M, 0.0E0, 1.0E0, U1, LDU1 )
         ENDIF
      ELSE
         CALL SLASCL( 'G', -1, -1, 1.0E0, W, M, N, A, LDA, INFO )
         IF ( INFO.NE.0 ) THEN
            RETURN
         END IF
         CALL SLACPY( 'A', M, N, A, LDA, WORK( IG11 ), LDG )
      END IF
*
      IF( PREPROCESSB ) THEN
         IWORK( 1:N ) = 0
         CALL SGEQP3( P, N, B, LDB, IWORK, WORK( ITAUB ),
     $                WORK( ITAUG ), LWORK - ITAUG + 1, INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
*        Determine rank
         ROWSB = 0
         DO I = 1, MIN( P, N )
            IF( ABS( B( I, I ) ).LE.ABSTOLB ) THEN
               EXIT
            END IF
            ROWSB = ROWSB + 1
         END DO
*        Copy full rank part into G
         CALL SLASET( 'L', ROWSB, N, 0.0E0, 0.0E0, WORK( IG21 ), LDG )
         CALL SLACPY( 'U', ROWSB, N, B, LDB, WORK( IG21 ), LDG )
         CALL SLAPMT( .FALSE., ROWSB, N, WORK( IG21 ), LDG, IWORK )
*        Initialize U2 although xORCSDB2BY1 will partially overwrite this
         IF( WANTU2 ) THEN
            CALL SLASET( 'A', P, P, 0.0E0, 1.0E0, U2, LDU2 )
         ENDIF
      ELSE
         CALL SLACPY( 'A', P, N, B, LDB, WORK( IG21 ), LDG )
      END IF
*
*     Compute the QR factorization with column pivoting GΠ = Q1 R1
*
      IWORK( 1:N ) = 0
      CALL SGEQP3( ROWSA + ROWSB, N, WORK( IG ), LDG, IWORK,
     $             WORK( ITAUG ),
     $             WORK( ISCRATCH ), LWORK - ISCRATCH + 1, INFO )
      IF( INFO.NE.0 ) THEN
         RETURN
      END IF
*
*     Compute the Frobenius norm of matrix G
*
      IF( NORMB.LE.UNFL ) THEN
         NORMG = W * NORMA
      ELSE
         NORMG = NORMB * SQRT( 1.0E0 + ( ( W * NORMA ) / NORMB )**2 )
      ENDIF
      ABSTOLG = TOL * MAX( ROWSA + ROWSB, N ) * MAX( NORMG, UNFL )
*
      IF( ISNAN(NORMG) ) THEN
         INFO = 103
         CALL XERBLA( 'SGGQRCS', INFO )
         RETURN
      ENDIF
*
*     Determine the rank of G
*
      RANK = 0
      DO I = 0, RANKMAXG - 1
         IF( ABS( WORK( IG + ( I * LDG + I ) ) ).LE.ABSTOLG ) THEN
            EXIT
         END IF
         RANK = RANK + 1
      END DO
*
*     Handle rank=0 case
*
      IF( RANK.EQ.0 ) THEN
         IF( WANTU1 ) THEN
            CALL SLASET( 'A', M, M, 0.0E0, 1.0E0, U1, LDU1 )
         END IF
         IF( WANTU2 ) THEN
            CALL SLASET( 'A', P, P, 0.0E0, 1.0E0, U2, LDU2 )
         END IF
*
         WORK( 1 ) = REAL( LWKOPT )
         RETURN
      END IF
*
*     Copy R1( 1:RANK, : ) into A, B
*
      IF( WANTX ) THEN
         IF( RANK.LE.M ) THEN
            CALL SLACPY( 'U', RANK, N, WORK( IG ), LDG, A, LDA )
         ELSE
            CALL SLACPY( 'U', M, N, WORK( IG ), LDG, A, LDA )
            CALL SLACPY( 'U', RANK - M, N - M, WORK( IG22 ), LDG,
     $                   B, LDB )
         END IF
      END IF
*
*     Explicitly form Q1 so that we can compute the CS decomposition
*
      CALL SORGQR( ROWSA + ROWSB, RANK, RANK, WORK( IG ), LDG,
     $             WORK( ITAUG ),
     $             WORK( ISCRATCH ), LWORK - ISCRATCH + 1, INFO )
      IF ( INFO.NE.0 ) THEN
         RETURN
      END IF
*
*     Compute the CS decomposition of Q1( :, 1:RANK )
*
      CALL SORCSD2BY1( JOBU1, JOBU2, JOBX, ROWSA + ROWSB, ROWSA, RANK,
     $                 WORK( IG11 ), LDG, WORK( IG21 ), LDG,
     $                 BETA,
     $                 U1, LDU1, U2, LDU2, X, LDX,
     $                 WORK( ITAUG ), LWORK - ITAUG + 1,
     $                 IWORK( N + 1 ), INFO )
      IF( INFO.NE.0 ) THEN
         RETURN
      END IF
*
*     Apply orthogonal factors of QR decomposition of A, B to U1, U2
*
      IF( PREPROCESSA .AND. WANTU1 ) THEN
         CALL SORMQR( 'RANK', 'N', M, M, ROWSA, A, LDA,
     $                WORK( ITAUA ), U1, LDU1,
     $                WORK( IG ), LWORK - IG + 1, INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         ENDIF
      ENDIF
*
      IF( PREPROCESSB .AND. WANTU2 ) THEN
         CALL SORMQR( 'RANK', 'N', P, P, ROWSB, B, LDB,
     $                WORK( ITAUB ), U2, LDU2,
     $                WORK( IG ), LWORK - IG + 1, INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         ENDIF
      ENDIF
*
*     Debug
*
      WORK( :LWORK ) = NAN
*
*     Compute X = V^T R1( 1:RANK, : ) and adjust for matrix scaling
*
      IF( WANTX ) THEN
         IF ( RANK.LE.M ) THEN
            CALL SGEMM( 'N', 'N', RANK, N - RANK, RANK,
     $                  1.0E0, X, LDX, A(1, RANK + 1), LDA,
     $                  0.0E0, X( 1, RANK + 1 ), LDX )
            CALL STRMM( 'R', 'U', 'N', 'N', RANK, RANK, 1.0E0,
     $                  A, LDA, X, LDX )
         ELSE
            CALL SLACPY( 'U', M, N, A, LDA, WORK( IG ), LDG )
            CALL SLACPY( 'U', RANK - M, N - M, B, LDB,
     $                   WORK( IG22 ), LDG )
            CALL SGEMM( 'N', 'N', RANK, N - RANK, RANK,
     $                  1.0E0, X, LDX,
     $                  WORK( IG + RANK * LDG ), LDG,
     $                  0.0E0, X( 1, RANK + 1 ), LDX )
            CALL STRMM( 'R', 'U', 'N', 'N', RANK, RANK, 1.0E0,
     $                  WORK( IG ), LDG, X, LDX )
         END IF
*        Revert column permutation Π by permuting the columns of X
         CALL SLAPMT( .FALSE., RANK, N, X, LDX, IWORK )
      END IF
*
*     Fix column order of U2
*     Because of the QR decomposition in the pre-processing, the first
*     rank(B) columns of U2 are a basis of range(B) but for matrix B,
*     the CS values are in ascending order. If B is singular, then the
*     first P - rank(B) columns should be a basis for the complement of
*     range(B). For this reason, the columns must be re-ordered.
*
      IF( PREPROCESSB .AND. WANTU2 .AND. ROWSB.LT.P ) THEN
         DO I = 1, ROWSB
            IWORK( I ) = P - ROWSB + I
         ENDDO
         IF( ROWSB.LE.P-ROWSB ) THEN
            DO I = ROWSB + 1, P - ROWSB
               IWORK( I ) = I
            ENDDO
            DO I = P - ROWSB + 1, P
               IWORK( I ) = I - (P - ROWSB)
            ENDDO
         ELSE
            DO I = ROWSB + 1, P
               IWORK( I ) = I - ROWSB
            ENDDO
         ENDIF
         CALL SLAPMT( .FALSE., P, P, U2, LDU2, IWORK )
      ENDIF
*
*     The pre-processing may reduce the number of angles that need to
*     be computed by the CS decomposition. Add these angles to the list
*     of angles.
*
      K = MIN( M, P, RANK, M + P - RANK )
      K1 = MAX( RANK - P, 0 )
      K2 = MAX( RANK - M, 0 )
*     Keep in mind the pre-processing might be disabled before
*     "optimizing" the expressions below
      KP = MIN( ROWSA, ROWSB, RANK, ROWSA + ROWSB - RANK )
      K1P = MAX( RANK - ROWSB, 0 )
      K2P = MAX( RANK - ROWSA, 0 )
*      PRINT*, "K , K1 , K2 ", K, K1, K2
*      PRINT*, "K', K1', K2'", KP, K1P, K2P
*     assert!(k == kp + k1p + k2p);
      IF( K.NE.KP + K1P - K1 + K2P - K2 ) THEN
        PRINT*, "k != k' + k1' - k1 + k2' - k2 !"
        PRINT*, "K, K1, K2   ", K, K1, K2
        PRINT*, "K', K1', K2'", KP, K1P, K2P
        INFO = 100
      ENDIF
*
      IF( K1P.GT.K1 ) THEN
*        Copy backwards because of a possible overlap
         DO I = KP, 1, -1
            BETA( I + K1P - K1 ) = BETA( I )
         ENDDO
         BETA( 1:K1P-K1 ) = 0.0E0
      ENDIF
      IF( K2P.GT.K2 ) THEN
         DO I = 1, K2P - K2
            BETA( I + K1P - K1 + KP ) = ACOS( 0.0E0 )
         ENDDO
      ENDIF
*
*     Adjust generalized singular values for matrix scaling
*     Compute sine, cosine values
*     Prepare row scaling of X
*
      DO I = 1, K
         THETA = BETA( I )
*        Do not adjust singular value if THETA is greater
*        than pi/2 (infinite singular values won't change)
         IF( COS( THETA ).LE.0.0E0 ) THEN
            ALPHA( I ) = 0.0E0
            BETA( I ) = 1.0E0
            IF( WANTX ) THEN
               WORK( I ) = 1.0E0
            END IF
         ELSE
*           iota comes in the greek alphabet after theta
            IOTA = ATAN( W * TAN( THETA ) )
*           ensure sine, cosine divisor is far away from zero
*           w is a power of two and will cause no trouble
            IF( SIN( IOTA ) .GE. COS( IOTA ) ) THEN
               ALPHA( I ) = ( SIN( IOTA ) / TAN( THETA ) ) / W
               BETA( I ) = SIN( IOTA )
               IF( WANTX ) THEN
                  WORK( I ) = SIN( THETA ) / SIN( IOTA )
               END IF
            ELSE
               ALPHA( I ) = COS( IOTA )
               BETA( I ) = SIN( IOTA )
               IF( WANTX ) THEN
                  WORK( I ) = COS( THETA ) / COS( IOTA ) / W
               END IF
            END IF
         END IF
      END DO
*     Adjust rows of X for matrix scaling
      IF( WANTX ) THEN
         DO J = 1, N
            DO I = 1, K1
               X( I, J ) = X( I, J ) / W
            END DO
            DO I = 1, K
               X( I + K1, J ) = X( I + K1, J ) * WORK( I )
            END DO
         END DO
      END IF
*
      WORK( 1 ) = REAL( LWKOPT )
      RETURN
*
*     End of SGGQRCS
*
      END
