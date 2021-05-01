* Copyright 2021 Christoph Conrads
      SUBROUTINE SGETRP( M, N, X, LDX, Y, LDY, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --
*
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            M, N, LDX, LDY, INFO
*     ..
*     .. Array Arguments ..
      REAL               X( LDX, * ), Y( LDY, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            B
*     64 byte cache line size divided by four bytes for a float
      PARAMETER          ( B = 16 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, K, L
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDX.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LDY.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ENDIF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGETRP', -INFO )
         RETURN
      END IF
*
*     Transpose matrix
*
      DO J = 1, N + B - 1, B
         DO I = 1, M + B - 1, B
            DO L = J, MIN( J + B, N )
               DO K = I, MIN( I + B, M )
                  Y( L, K ) = X( K, L )
               ENDDO
            ENDDO
         ENDDO
      ENDDO
*
      RETURN
*
*     End of SGETRP
*
      END
