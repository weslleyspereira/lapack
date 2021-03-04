!> \brief \b SLANGE returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL(wp)         FUNCTION SLANGE( NORM, M, N, A, LDA, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM
!       INTEGER            LDA, M, N
!       ..
!       .. Array Arguments ..
!       REAL(wp)           A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLANGE  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> real matrix A.
!> \endverbatim
!>
!> \return SLANGE
!> \verbatim
!>
!>    SLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!>             (
!>             ( norm1(A),         NORM = '1', 'O' or 'o'
!>             (
!>             ( normI(A),         NORM = 'I' or 'i'
!>             (
!>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!>
!> where  norm1  denotes the  one norm of a matrix (maximum column sum),
!> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!> normF  denotes the  Frobenius norm of a matrix (square root of sum of
!> squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NORM
!> \verbatim
!>          NORM is CHARACTER*1
!>          Specifies the value to be returned in SLANGE as described
!>          above.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.  When M = 0,
!>          SLANGE is set to zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.  When N = 0,
!>          SLANGE is set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL(wp) array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(M,1).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL(wp) array, dimension (MAX(1,LWORK)),
!>          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
!>          referenced.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date December 2016
!
!> \ingroup realGEauxiliary
!
!> \par Contributors:
!  ==================
!>
!> Edward Anderson, Lockheed Martin
!> Weslley Pereira, University of Colorado Denver, USA
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Anderson E. (2017)
!>  Algorithm 978: Safe Scaling in the Level 1 BLAS
!>  ACM Trans Math Softw 44:1--28
!>  https://doi.org/10.1145/3061665
!>
!> \endverbatim
!
!  =====================================================================
function SLANGE( norm, m, n, a, lda, work )
   use LA_CONSTANTS, only: wp=>sp, zero=>szero, one=>sone
   use LA_XMAXNAN
   real(wp) :: SLANGE
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!  .. Scalar Arguments ..
   character :: norm
   integer :: lda, m, n
!  ..
!  .. Array Arguments ..
   real(wp) :: a( lda, * ), work( * )
!  ..
!  .. Local Scalars ..
   integer :: i, j
   real(wp) :: scale, sum, value
!  ..
!  .. External Functions ..
   logical :: LSAME
   external :: LSAME
!  ..
!  .. External Subroutines ..
   external :: SLASSQ
!  ..
!  .. Executable Statements ..
!
   if( min( m, n ) == 0 ) then
      value = zero
   else if( LSAME( norm, 'M' ) ) then
!
!     Find max(abs(A(i,j))).
!
      value = zero
      do j = 1, n
         do i = 1, m
            value = LA_MAXNAN( value, abs( a( i, j ) ) )
         end do
      end do
   else if( LSAME( norm, 'O' ) .or. norm == '1' ) then
!
!     Find norm1(A).
!
      value = zero
      do j = 1, n
         sum = zero
         do i = 1, m
            sum = sum + abs( a( i, j ) )
         end do
         value = LA_MAXNAN( value, sum )
      end do
   else if( LSAME( norm, 'I' ) ) then
!
!     Find normI(A).
!
      do i = 1, m
         work( i ) = zero
      end do
      do j = 1, n
         do i = 1, m
            work( i ) = work( i ) + abs( a( i, j ) )
         end do
      end do
      value = zero
      do i = 1, m
         value = LA_MAXNAN( value, work( i ) )
      end do
   else if( LSAME( norm, 'F' ) .or. LSAME( norm, 'E' ) ) then
!
!     Find normF(A).
!
      scale = one
      sum = zero
      do j = 1, n
         call SLASSQ( m, a( 1, j ), 1, scale, sum )
      end do
      value = scale*sqrt( sum )
   end if
!
   SLANGE = value
   return
end function SLANGE
