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
*       SUBROUTINE SGGQRCS( JOBU1, JOBU2, JOBX,
*      $                    HINTPREPA, HINTPREPB,
*      $                    M, N, P, RANK,
*      $                    SWAPPED,
*      $                    A, LDA, B, LDB,
*      $                    ALPHA, BETA,
*      $                    U1, LDU1, U2, LDU2, X, LDX,
*      $                    TOL,
*      $                    WORK, LWORK,
*      $                    IWORK, INFO )
*
*       .. Scalar Arguments ..
*       LOGICAL            SWAPPED
*       CHARACTER          JOBU1, JOBU2, JOBX, HINTPREPA, HINTPREPB
*       INTEGER            INFO, LDA, LDB, LDU1, LDU2, LDX,
*      $                   M, N, P, RANK, LWORK
*       REAL               TOL
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       REAL               A( LDA, * ), B( LDB, * ),
*      $                   ALPHA( N ), BETA( N ),
*      $                   U1( LDU1, * ), U2( LDU2, * ), X( LDX, * )
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
*> \param[in,out] HINTPREPA
*> \verbatim
*>          HINTPREPA is CHARACTER*1
*>          = 'Y':  Advises SGGQRCS to pre-process A;
*>          = 'N':  advises SGGQRCS not to pre-process A;
*>          = '?':  pre-processing decision left to algorithm.
*>
*>          On exit, the variable will be set to 'Y' if the matrix
*>          was pre-processed, 'N' otherwise.
*>          The variable may be modified by workspace queries.
*>
*>          See below for a discussion of pre-processing.
*> \endverbatim
*>
*> \param[in,out] HINTPREPB
*> \verbatim
*>          HINTPREPB is CHARACTER*1
*>          = 'Y':  Advises SGGQRCS to pre-process B;
*>          = 'N':  advises SGGQRCS not to pre-process B;
*>          = '?':  pre-processing decision left to algorithm.
*>
*>          On exit, the variable will be set to 'Y' if the matrix
*>          was pre-processed, 'N' otherwise.
*>          The variable may be modified by workspace queries.
*>
*>          See below for a discussion of pre-processing.
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
*>          SWAPPED is LOGICAL
*>          On exit, SWAPPED is true if SGGQRCS swapped the input
*>          matrices A, B and computed the GSVD of (B, A); false
*>          otherwise.
*>
*>          Swapping the matrices internally is necessary to achieve
*>          uncoditional backward stability.
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
*> \param[out] X
*> \verbatim
*>          X is REAL array, dimension (LDX,N)
*>          If JOBX = 'Y', X contains the RANKxN matrix X.
*>          If JOBX = 'N', X is not referenced.
*> \endverbatim
*>
*> \param[in] LDX
*> \verbatim
*>          LDX is INTEGER
*>          The leading dimension of the array X. LDX >= max(1,RANKMAXG)
*>          if JOBX = 'Y'; LDX >= 1 otherwise. RANKMAXG is the largest
*>          possible rank of the matrix G = (A**T, W*B**T)**T:
*>            RANKMAXG = MIN( MIN( M, N ) + MIN( P, N ), N ).
*> \endverbatim
*>
*> \param[in,out] TOL
*> \verbatim
*>          TOL is REAL
*>          This user-provided tolerance is used for the rank determination
*>          of the matrix
*>          * A if pre-processing of A is enabled,
*>          * B if pre-processing of B is enabled,
*>          * G = (A**T, W*B**T)**T.
*>          See the documentation of ABSTOL for details.
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
*>          > 0:  One of the subroutines failed. Its name will be
*>                printed by XERBLA.
*>          101:  The norm of A is not a number of infinite.
*>          102:  The norm of B is not a number of infinite.
*>          103:  The norm of G is not a number of infinite.
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
*>  ABSTOL
*>  ABSTOLA REAL
*>  ABSTOLB REAL
*>  ABSTOLG REAL
*>          ABSTOL are threshold values to determine the numerical rank
*>          of matrices. For a matrix Z, it is set to
*>              ABSTOLZ = TOL * MAX( #ROWS, #COLUMNS ) * norm(Z),
*>          where norm(Z) is the Frobenius norm of Z.
*>          The size of ABSTOL may affect the size of backward error of the
*>          decomposition.
*>
*>  ITAUGL  INTEGER
*>          The index of the first scalar factor of the elementary
*>          reflectors of the LQ decomposition of G.
*>
*>  ITAUGR  INTEGER
*>          The index of the first scalar factor of the elementary
*>          reflectors of the QR decomposition of G.
*>
*>  MAYBEPREPG LOGICAL
*>          Pre-processing the matrix G requires dedicated workspace.
*>          This variable captures the need to reserve this workspace
*>          but it is not guaranteed that G will be pre-processed.
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
*>  SGGQRCS can pre-process the input matrices A and B. Pre-processing
*>  can noticeably sped up the GSVD computation and it limits the
*>  workspace size for matrices with more rows than columns but it does
*>  not influence the accuracy of the results. Hence, the pre-processing
*>  is best thought of as a possibility to optimize performance. If the
*>  matrix rank is small in comparison to the number of its rows, then
*>  it is recommended to pre-process the matrix and to pass 'Y' as a
*>  pre-processing hint. If the matrix has (almost) full row rank,
*>  then the pre-processing may not speed up the GSVD computation and it
*>  is best to pass 'N' as a pre-processing hint. In cases where the
*>  matrix rank is not known or if the user does not care, '?' should be
*>  passed as hint. SGGQRCS may ignore these _hints_, in particular for
*>  matrices with much more rows than columns the pre-processing is
*>  necessary to bound the workspace size by a multiple of the number of
*>  columns. This design choice allows LAPACK implementations to
*>  customize their pre-processing criteria.
*>
*>  SGGQRCS should be significantly faster than SGGSVD3 for large
*>  matrices because the matrices A and B are reduced to a pair of
*>  well-conditioned bidiagonal matrices instead of pairs of upper
*>  triangular matrices. On the downside, SGGQRCS requires a much larger
*>  workspace whose dimension must be queried at run-time. SGGQRCS also
*>  offers no guarantees which of the two possible diagonal matrices
*>  is used for the matrix factorization.
*>
*> \par Notes to Implementors
*  ==========================
*>
*>  * The workspace queries of xORGLQ and xORMLQ rely on the matrix G
*>    having at least twice as many columns as rows.
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
      LOGICAL            MAYBEPREPG,
     $                   PREPROCESSA, PREPROCESSB, PREPROCESSG,
     $                   WANTU1, WANTU2, WANTX, LQUERY
      INTEGER            I, J,
     $                   K, K1, K2, KP, K1P, K2P,
     $                   COLS,
     $                   RANKMAXA, RANKMAXB, RANKMAXG, ROWSA, ROWSB,
     $                   ITAUA, ITAUB, ITAUGL, ITAUGR, IMAT, ISCRATCH,
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
      EXTERNAL           SGEMM, SGEQP3, SGETRP, SLACPY, SLAPMT, SLASCL,
     $                   SLASET, SORGQR, SORCSD2BY1, SORMLQ, SORMQR,
     $                   XERBLA
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
      WANTU1 = LSAME( JOBU1, 'Y' )
      WANTU2 = LSAME( JOBU2, 'Y' )
      WANTX = LSAME( JOBX, 'Y' )
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
      PREPROCESSG = 2 * ( M + P ).LE.N
*
      RANK = -1
      COLS = N
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
      IF( PREPROCESSG ) THEN
         LDG = MAX( N, 1 )
      ELSE
         LDG = MAX( ROWSA + ROWSB, 1 )
      ENDIF
*     Compute offsets into workspace
      ITAUA = 1
      IF( PREPROCESSA ) THEN
         ITAUB = ITAUA + RANKMAXA
      ELSE
         ITAUB = ITAUA
      ENDIF
      IF( PREPROCESSB ) THEN
         IMAT = ITAUB + RANKMAXB
      ELSE
         IMAT = ITAUB
      ENDIF
      ITAUGL = IMAT + MAX( ROWSA + ROWSB, 1 ) * N
*     TODO: reserving workspace for the scalar factors of the LQ
*     factorization with row pivoting is only needed if WANTX
      IF( PREPROCESSG ) THEN
         ITAUGR = ITAUGL + ROWSA + ROWSB
      ELSE
         ITAUGR = ITAUGL
      ENDIF
      ISCRATCH = ITAUGR + RANKMAXG
      IG11 = IMAT
      IG21 = IMAT + ROWSA
      IG22 = IMAT + LDG * M + M
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
         CALL SGGQRCS( JOBU2, JOBU1, JOBX, HINTPREPB, HINTPREPA,
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
      IF( PREPROCESSG ) THEN
         CALL SGEQP3( N, M + P, WORK( IMAT ), LDG, IWORK, WORK,
     $                WORK, -1, INFO )
         LWKMIN = MAX( LWKMIN, 3 * ( M + P ) + 1 + ITAUB )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) )  + ITAUB )
      ENDIF
*
      IF( PREPROCESSA ) THEN
         CALL SGEQP3( M, N, A, LDA, IWORK, WORK( ITAUA ),
     $                WORK, -1, INFO )
         LWKMIN = MAX( LWKMIN, 3 * N + 1        + ITAUB )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) + ITAUB )
      ENDIF
*
      IF( PREPROCESSB ) THEN
         CALL SGEQP3( P, N, B, LDB, IWORK, WORK( ITAUB ),
     $                WORK, -1, INFO )
         LWKMIN = MAX( LWKMIN, 3 * N + 1        + IMAT )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) + IMAT )
      ENDIF
*
      IF( PREPROCESSG ) THEN
         CALL SGEQRF( ROWSA + ROWSB, RANKMAXG, WORK, LDG, WORK,
     $                WORK, -1, INFO )
         LWKMIN = MAX( LWKMIN, INT( WORK( 1 ) ) + ISCRATCH )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) + ISCRATCH )
      ELSE
         CALL SGEQP3( ROWSA + ROWSB, N, WORK( IMAT ), LDG, IWORK,
     $             WORK( ITAUGR ), WORK, -1, INFO )
         LWKMIN = MAX( LWKMIN, 3 * N + 1        + ISCRATCH )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) + ISCRATCH )
      ENDIF
*
      CALL SORGQR( ROWSA + ROWSB, RANKMAXG, RANKMAXG, WORK( IMAT ), LDG,
     $             WORK( ITAUGR ), WORK, -1, INFO )
      LWKMIN = MAX( LWKMIN, RANKMAXG         + ISCRATCH )
      LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) + ISCRATCH )
*
      CALL SORCSD2BY1( JOBU1, JOBU2, JOBX,
     $                 ROWSA + ROWSB, ROWSA, RANKMAXG,
     $                 WORK( IG11 ), LDG, WORK( IG21 ), LDG,
     $                 BETA,
     $                 U1, LDU1, U2, LDU2, X, LDX,
     $                 WORK, -1, IWORK, INFO )
      LWKMIN = MAX( LWKMIN, INT( WORK( 1 ) ) + ITAUGR )
      LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) + ITAUGR )
*
      IF( PREPROCESSG .AND. WANTX ) THEN
         I = ITAUGR - ITAUGL
         CALL SORMQR( 'R', 'T', I, N, I, WORK, LDG, WORK, X, LDX,
     $                WORK, -1, INFO )
         LWKMIN = MAX( LWKMIN, MAX( 1, N )      + ITAUGR )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) + ITAUGR )
      ENDIF
*
      IF( PREPROCESSA .AND. WANTU1 ) THEN
         CALL SORMQR( 'L', 'N', M, M, RANKMAXA, A, LDA,
     $                WORK, U1, LDU1,
     $                WORK, -1, INFO )
         LWKMIN = MAX( LWKMIN, MAX( 1, N )      + IMAT )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) + IMAT )
      ENDIF
*
      IF( PREPROCESSB .AND. WANTU2 ) THEN
         CALL SORMQR( 'L', 'N', P, P, RANKMAXB, B, LDB,
     $                WORK, U2, LDU2,
     $                WORK, -1, INFO )
         LWKMIN = MAX( LWKMIN, MAX( 1, N )      + IMAT )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) + IMAT )
      ENDIF
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
*     Compute the Frobenius norm of matrix G
*
      IF( NORMB.LE.UNFL ) THEN
         NORMG = W * NORMA
      ELSE
         NORMG = NORMB * SQRT( 1.0E0 + ( ( W * NORMA ) / NORMB )**2 )
      ENDIF
*
      ABSTOLG = TOL * MAX( ROWSA + ROWSB, N ) * MAX( NORMG, UNFL )
*
      IF( ISNAN(NORMG) ) THEN
         INFO = 103
         CALL XERBLA( 'SGGQRCS', INFO )
         RETURN
      ENDIF
*
      CALL SLASCL( 'G', -1, -1, 1.0E0, W, M, N, A, LDA, INFO )
      IF ( INFO.NE.0 ) THEN
         RETURN
      ENDIF
*
*     Remove unnecessary matrix columns
*
      IF( PREPROCESSG ) THEN
*        assert!( 2 * (m + p) <= n );
         CALL SGETRP( M, N, A, LDA, WORK( IMAT ), LDG, INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         ENDIF
         CALL SGETRP( P, N, B, LDB, WORK( IMAT + M * LDG ), LDG, INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         ENDIF
*
         IWORK( 1:M+P ) = 0
         CALL SGEQP3( N, M + P, WORK( IMAT ), LDG,
     $                IWORK, WORK( ITAUGL ),
     $                WORK( ITAUGR ), LWORK - ITAUGR + 1, INFO )
*        Determine the rank of G
         RANK = 0
         DO I = 0, M + P - 1
            IF( ABS( WORK( IMAT + ( I * LDG + I ) ) ).LE.ABSTOLG ) THEN
               EXIT
            END IF
            RANK = RANK + 1
         END DO
         COLS = RANK
*        Multiply orthogonal factor with A
         CALL SORMQR( 'R', 'N', M, N, RANK, WORK( IMAT ), LDG,
     $                WORK( ITAUGL ), A, LDA,
     $                WORK( ITAUGR ), LWORK - ITAUGR + 1, INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         ENDIF
*        Multiply orthogonal factor with B
         CALL SORMQR( 'R', 'N', P, N, RANK, WORK( IMAT ), LDG,
     $                WORK( ITAUGL ), B, LDB,
     $                WORK( ITAUGR ), LWORK - ITAUGR + 1, INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         ENDIF
*        Save upper half of elementary reflectors; this half will be
*        re-used  by G
         IF( WANTX ) THEN
            CALL SLACPY( 'L', ROWSA + ROWSB, COLS, WORK( IMAT ), LDG,
     $                   X( 1, RANK + 1 ), LDX )
         ENDIF
      ENDIF
*
*     Attempt to remove unnecessary matrix rows
*     Copy matrices A, B or their full-rank factors, respectively, into
*     the LDG x N matrix G
*
      IF( PREPROCESSA ) THEN
         IWORK( 1:COLS ) = 0
         CALL SGEQP3( M, COLS, A, LDA, IWORK, WORK( ITAUA ),
     $                WORK( ITAUGR ), LWORK - ITAUGR + 1, INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
*        Determine rank of A
         ROWSA = 0
         DO I = 1, MIN( M, COLS )
            IF( ABS( A( I, I ) ).LE.ABSTOLA ) THEN
               EXIT
            ENDIF
            ROWSA = ROWSA + 1
         END DO
         IG21 = IG11 + ROWSA
*        Copy full rank part into G
         CALL SLASET( 'L', ROWSA, COLS, 0.0E0, 0.0E0,
     $                WORK( IG11 ), LDG )
         CALL SLACPY( 'U', ROWSA, COLS, A, LDA, WORK( IG11 ), LDG )
         CALL SLAPMT( .FALSE., ROWSA, COLS, WORK( IG11 ), LDG, IWORK )
*        Initialize U1 although xORCSDB2BY1 will partially overwrite this
         IF( WANTU1 ) THEN
             CALL SLASET( 'A', M, M, 0.0E0, 1.0E0, U1, LDU1 )
         ENDIF
      ELSE
         CALL SLACPY( 'A', M, COLS, A, LDA, WORK( IG11 ), LDG )
      END IF
*
      IF( PREPROCESSB ) THEN
         IWORK( 1:COLS ) = 0
         CALL SGEQP3( P, COLS, B, LDB, IWORK, WORK( ITAUB ),
     $                WORK( ITAUGR ), LWORK - ITAUGR + 1, INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
*        Determine rank
         ROWSB = 0
         DO I = 1, MIN( P, COLS )
            IF( ABS( B( I, I ) ).LE.ABSTOLB ) THEN
               EXIT
            END IF
            ROWSB = ROWSB + 1
         END DO
*        Copy full rank part into G
         CALL SLASET( 'L', ROWSB, COLS, 0.0E0, 0.0E0,
     $                WORK( IG21 ), LDG )
         CALL SLACPY( 'U', ROWSB, COLS, B, LDB, WORK( IG21 ), LDG )
         CALL SLAPMT( .FALSE., ROWSB, COLS, WORK( IG21 ), LDG, IWORK )
*        Initialize U2 although xORCSDB2BY1 will partially overwrite this
         IF( WANTU2 ) THEN
            CALL SLASET( 'A', P, P, 0.0E0, 1.0E0, U2, LDU2 )
         ENDIF
      ELSE
         CALL SLACPY( 'A', P, COLS, B, LDB, WORK( IG21 ), LDG )
      END IF
*
*     Compute QR decomposition of G
*
      IF( PREPROCESSG ) THEN
*        assert!( COLS < N );
         CALL SGEQRF( ROWSA + ROWSB, COLS, WORK( IG11 ), LDG,
     $                WORK( ITAUGR ),
     $                WORK( ISCRATCH ), LWORK - ISCRATCH + 1, INFO )
      ELSE
*        assert!( COLS == N );
*        assert!( RANK < 0 );
*
*        Compute the QR factorization with column pivoting GΠ = Q1 R1
         IWORK( 1:N ) = 0
         CALL SGEQP3( ROWSA + ROWSB, COLS, WORK( IG11 ), LDG, IWORK,
     $                WORK( ITAUGR ),
     $                WORK( ISCRATCH ), LWORK - ISCRATCH + 1, INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
*        Determine the rank of G
         RANK = 0
         DO I = 0, MIN( ROWSA + ROWSB, N ) - 1
            IF( ABS( WORK( IMAT + ( I * LDG + I ) ) ).LE.ABSTOLG ) THEN
               EXIT
            END IF
            RANK = RANK + 1
         END DO
      ENDIF
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
            CALL SLACPY( 'U', RANK, COLS, WORK( IG11 ), LDG, A, LDA )
         ELSE
            CALL SLACPY( 'U', M, COLS, WORK( IG11 ), LDG, A, LDA )
            CALL SLACPY( 'U', RANK - M, COLS - M, WORK( IG22 ), LDG,
     $                   B, LDB )
         END IF
      END IF
*     DEBUG
      IF( PREPROCESSG ) THEN
          CALL SLASET( 'U', N, COLS, NAN, NAN, WORK( IG11 ), LDG )
      ELSE
          CALL SLASET( 'U', RANK, COLS, NAN, NAN, WORK( IG11 ), LDG )
      ENDIF
*
*     Explicitly form Q1 so that we can compute the CS decomposition
*
*
      CALL SORGQR( ROWSA + ROWSB, RANK, RANK, WORK( IG11 ), LDG,
     $             WORK( ITAUGR ),
     $             WORK( ISCRATCH ), LWORK - ISCRATCH + 1, INFO )
      IF ( INFO.NE.0 ) THEN
         RETURN
      END IF
*     DEBUG
      WORK( ITAUGR:LWORK ) = NAN
*
*     The pre-processing may reduce the number of angles that need to
*     be computed by the CS decomposition:
*     * angles 1, 2, ..., K1P-K1 are known to be zero,
*     * angles K1P-K1+1, ..., K1P-K1+KP are computed by xORCSD2BY1,
*     * angles K1P-K1+KP+1, ..., K are known to be pi/2.
*
*     It holds that
*     * K = KP + (K1P - K1) + (K2P - K2),
*     * RANK = K + K1 + K2 = KP + K1P + K2P
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
*      PRINT*, "MNP         ", RANK, K + K1 + K2, KP + K1P + K2P
*     assert!(k == kp + k1p + k2p);
      IF( K.NE.KP + K1P - K1 + K2P - K2 ) THEN
        PRINT*, "k != k' + k1' - k1 + k2' - k2 !"
        PRINT*, "K, K1, K2   ", K, K1, K2
        PRINT*, "K', K1', K2'", KP, K1P, K2P
        INFO = 100
        RETURN
      ENDIF
      IF( RANK.NE.KP + K1P + K2P) THEN
        INFO = 99
        CALL XERBLA( 'SGGQRCS', INFO )
        RETURN
      ENDIF
*
*     Compute the CS decomposition of Q1( :, 1:RANK )
*
      CALL SORCSD2BY1( JOBU1, JOBU2, JOBX, ROWSA + ROWSB, ROWSA, RANK,
     $                 WORK( IG11 ), LDG, WORK( IG21 ), LDG,
     $                 BETA( K1P - K1 + 1 ),
     $                 U1, LDU1, U2, LDU2, X, LDX,
     $                 WORK( ITAUGR ), LWORK - ITAUGR + 1,
     $                 IWORK( COLS + 1 ), INFO )
      IF( INFO.NE.0 ) THEN
         RETURN
      END IF
*     DEBUG
      WORK( ITAUGR:LWORK ) = NAN
*
*     Compute X = ( V^T R1( 1:RANK, : ) ) Q, where Q is the orthogonal
*     factor of the LQ factorization of G if pre-processing was enabled
*     and the identity matrix otherwise. V and R1 should be multiplied
*     first to minimize the FLOP count.
*
      IF( WANTX ) THEN
*        Move upper half of elementary reflectors of LQ back
         IF( PREPROCESSG ) THEN
            CALL SLACPY( 'L', ROWSA + ROWSB, COLS, X( 1,RANK + 1 ), LDX,
     $                   WORK( IG11 ), LDG )
*           Debug
            CALL SLASET( 'A', RANK, N - RANK, NAN, NAN,
     $                   X( 1, RANK + 1 ), LDX )
         ENDIF
*
         IF ( RANK.LE.M ) THEN
            IF( COLS.GT.RANK ) THEN
               CALL SGEMM( 'N', 'N', RANK, COLS - RANK, RANK,
     $                     1.0E0, X, LDX, A( 1, RANK + 1 ), LDA,
     $                     0.0E0, X( 1, RANK + 1 ), LDX )
            END IF
            CALL STRMM( 'R', 'U', 'N', 'N', RANK, RANK, 1.0E0,
     $                  A, LDA, X, LDX )
         ELSE
            CALL SLACPY( 'U', M, COLS, A, LDA, WORK( IG11 ), LDG )
            CALL SLACPY( 'U', RANK - M, COLS - M, B, LDB,
     $                   WORK( IG22 ), LDG )
            IF( COLS.GT.RANK ) THEN
               CALL SGEMM( 'N', 'N', RANK, COLS - RANK, RANK,
     $                     1.0E0, X, LDX,
     $                     WORK( IG11 + RANK * LDG ), LDG,
     $                     0.0E0, X( 1, RANK + 1 ), LDX )
            ENDIF
            CALL STRMM( 'R', 'U', 'N', 'N', RANK, RANK, 1.0E0,
     $                  WORK( IG11 ), LDG, X, LDX )
         END IF
*
         IF( PREPROCESSG ) THEN
            CALL SLASET( 'G', RANK, N - RANK, 0.0E0, 0.0E0,
     $                   X( 1, RANK + 1 ), LDX )
            CALL SORMQR( 'R', 'T', RANK, N, RANK, WORK( IMAT ), LDG,
     $                   WORK( ITAUGL ), X, LDX,
     $                   WORK( ITAUGR ), LWORK - ITAUGR + 1, INFO )
         ELSE
*           Revert column permutation Π by permuting the columns of X
            CALL SLAPMT( .FALSE., RANK, COLS, X, LDX, IWORK )
         ENDIF
      ENDIF
*     DEBUG
      WORK( IMAT:LWORK ) = NAN
*
*     Apply orthogonal factors of QR decomposition of A, B to U1, U2
*
      IF( PREPROCESSA .AND. WANTU1 ) THEN
         CALL SORMQR( 'L', 'N', M, M, ROWSA, A, LDA,
     $                WORK( ITAUA ), U1, LDU1,
     $                WORK( IMAT ), LWORK - IMAT + 1, INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         ENDIF
      ENDIF
*
      IF( PREPROCESSB .AND. WANTU2 ) THEN
         CALL SORMQR( 'L', 'N', P, P, ROWSB, B, LDB,
     $                WORK( ITAUB ), U2, LDU2,
     $                WORK( IMAT ), LWORK - IMAT + 1, INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         ENDIF
*        Fix column order of U2
*        Because of the QR decomposition in the pre-processing, the first
*        rank(B) columns of U2 are a basis of range(B) but for matrix B,
*        the CS values are in ascending order. If B does not have full
*        row rank, then the first P - rank(B) columns should be a basis
*        for the complement of range(B) and the columns of U2 must be
*        re-ordered.
*        Ensure the column permutation of G was reverted before
*        overwriting IWORK.
         IF( ROWSB.LT.P ) THEN
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
      ENDIF
*
*     Debug
*
      WORK( :LWORK ) = NAN
*
*     Adjust generalized singular values for matrix scaling
*     Compute sine, cosine values
*     Prepare row scaling of X
*
*     Reading hints:
*     * the first K1 angles are known to be zero by the caller. thus,
*       the sine, cosine, and row scaling information associated with
*       the I-th angle, I = 1, 2, ..., RANK,
*       * is not stored for I = 1, 2, ..., K1,
*       * at index I - K1 for I = K1+1, K1+2, ..., K,
*       * is not stored for K+1, K+2, ..., RANK.
*     * row scaling factors for X will not be stored if they won't be
*       used (e.g., because they are known to be one).
*     * the row scaling factors for X will be computed irrespective of
*       the value of WANTX. this improves readability.
*     * the code below this paragraph is one large iteration of I from 1
*       to K with different DO statements for the different cases
      DO I = 1, K1P - K1
         ALPHA( I ) = 1.0E0
         BETA( I ) = 0.0E0
*        DEBUG
         WORK( I ) = NAN
      ENDDO
      DO I = K1P - K1 + 1, K1P - K1 + KP
         THETA = BETA( I )
*        Do not adjust singular value if THETA is greater
*        than pi/2 (infinite singular values won't change)
         IF( COS( THETA ).LE.0.0E0 ) THEN
            ALPHA( I ) = 0.0E0
            BETA( I ) = 1.0E0
            WORK( I ) = 1.0E0
         ELSE
*           iota comes in the greek alphabet after theta
            IOTA = ATAN( W * TAN( THETA ) )
*           ensure sine, cosine divisor is far away from zero
*           w is a power of two and will cause no trouble
            IF( SIN( IOTA ) .GE. COS( IOTA ) ) THEN
               ALPHA( I ) = ( SIN( IOTA ) / TAN( THETA ) ) / W
               BETA( I ) = SIN( IOTA )
               WORK( I ) = SIN( THETA ) / SIN( IOTA )
            ELSE
               ALPHA( I ) = COS( IOTA )
               BETA( I ) = SIN( IOTA )
               WORK( I ) = COS( THETA ) / COS( IOTA ) / W
            END IF
         END IF
      END DO
      DO I = K1P - K1 + KP + 1, K
         ALPHA( I ) = 0.0E0
         BETA( I ) = 1.0E0
*        DEBUG
         WORK( I ) = NAN
      ENDDO
*     Adjust rows of X for matrix scaling
      IF( WANTX ) THEN
         DO J = 1, N
            DO I = 1, K1P
               X( I, J ) = X( I, J ) / W
            END DO
            DO I = K1P + 1, K1P + KP
               X( I, J ) = X( I, J ) * WORK( I - K1 )
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
