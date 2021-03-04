module LA_XMAXNAN
!
!  LAPACK auxiliary routine
!  E. Anderson
!  October 30, 2020
!
   interface LA_MAXNAN
      module procedure SMAXNAN
      module procedure DMAXNAN
   end interface
contains
   
function SMAXNAN( x, y )
   use LA_CONSTANTS, only: wp=>sp
   use LA_XISNAN
   real(wp) :: SMAXNAN
!
!  Purpose
!  =======
!
!  SMAXNAN is a NaN-aware version of the max function that
!  returns NaN if either x or y is NaN, otherwise it returns
!  max(x,y).
!
!  Arguments
!  =========
!
!  X       (input) REAL
!          The scalar x.
!
!  Y       (input) REAL
!          The scalar y.
!
! =====================================================================
!
!  .. Scalar Arguments ..
   real(wp) :: x, y
!  ..
   if( LA_ISNAN(x) ) then
      SMAXNAN = x
   else if( LA_ISNAN(y) ) then
      SMAXNAN = y
   else
      SMAXNAN = max( x, y )
   end if
end function SMAXNAN

function DMAXNAN( x, y )
   use LA_CONSTANTS, only: wp=>dp
   use LA_XISNAN
   real(wp) :: DMAXNAN
!
!  Purpose
!  =======
!
!  DMAXNAN is a NaN-aware version of the max function that
!  returns NaN if either x or y is NaN, otherwise it returns
!  max(x,y).
!
!  Arguments
!  =========
!
!  X       (input) REAL
!          The scalar x.
!
!  Y       (input) REAL
!          The scalar y.
!
! =====================================================================
!
!  .. Scalar Arguments ..
   real(wp) :: x, y
!  ..
   if( LA_ISNAN(x) ) then
      DMAXNAN = x
   else if( LA_ISNAN(y) ) then
      DMAXNAN = y
   else
      DMAXNAN = max( x, y )
   end if
end function DMAXNAN

end module LA_XMAXNAN
