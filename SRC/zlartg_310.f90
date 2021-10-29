!> \brief \b ZLARTG generates a plane rotation with real cosine and complex sine.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLARTG( F, G, C, S, R )
!
!       .. Scalar Arguments ..
!       REAL(wp)           C
!       COMPLEX(wp)        F, G, R, S
!       ..
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLARTG generates a plane rotation so that
!>
!>    [  C         S  ] . [ F ]  =  [ R ]
!>    [ -conjg(S)  C  ]   [ G ]     [ 0 ]
!>
!> where C is real and C**2 + |S|**2 = 1.
!>
!> The mathematical formulas used for C and S are
!>
!>    sgn(x) = {  x / |x|,   x != 0
!>             {  1,         x = 0
!>
!>    R = sgn(F) * sqrt(|F|**2 + |G|**2)
!>
!>    C = |F| / sqrt(|F|**2 + |G|**2)
!>
!>    S = sgn(F) * conjg(G) / sqrt(|F|**2 + |G|**2)
!>
!> When F and G are real, the formulas simplify to C = F/R and
!> S = G/R, and the returned values of C, S, and R should be
!> identical to those returned by DLARTG.
!>
!> The algorithm used to compute these quantities incorporates scaling
!> to avoid overflow or underflow in computing the square root of the
!> sum of squares.
!>
!> This is a faster version of the BLAS1 routine ZROTG, except for
!> the following differences:
!>    F and G are unchanged on return.
!>    If G=0, then C=1 and S=0.
!>    If F=0, then C=0 and S is chosen so that R is real.
!>
!> Below, wp=>dp stands for double precision from LA_CONSTANTS module.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] F
!> \verbatim
!>          F is COMPLEX(wp)
!>          The first component of vector to be rotated.
!> \endverbatim
!>
!> \param[in] G
!> \verbatim
!>          G is COMPLEX(wp)
!>          The second component of vector to be rotated.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is REAL(wp)
!>          The cosine of the rotation.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is COMPLEX(wp)
!>          The sine of the rotation.
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is COMPLEX(wp)
!>          The nonzero component of the rotated vector.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Edward Anderson, Lockheed Martin
!
!> \date August 2016
!
!> \ingroup OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
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
subroutine ZLARTG( f, g, c, s, r )
   use LA_CONSTANTS, &
   only: wp=>dp, zero=>dzero, one=>done, two=>dtwo, czero=>zzero, &
         rtmin=>drtmin, rtmax=>drtmax, safmin=>dsafmin, safmax=>dsafmax
!
!  -- LAPACK auxiliary routine (version 3.10.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     February 2021
!
!  .. Scalar Arguments ..
   real(wp)           c
   complex(wp)        f, g, r, s
!  ..
!  .. Local Scalars ..
   real(wp) :: d, f1, f2, g1, g2, h2, u, v, w, noise1, noise2
   complex(wp) :: fs, gs, t
!  ..
!  .. Intrinsic Functions ..
   intrinsic :: abs, aimag, conjg, max, min, real, sqrt
!  ..
!  .. Statement Functions ..
   real(wp) :: ABSSQ
!  ..
!  .. Statement Function definitions ..
   ABSSQ( t ) = real( t )**2 + aimag( t )**2
!  ..
!  .. Executable Statements ..
!
   if( g == czero ) then
      stop 1
      c = one
      s = czero
      r = f
   else if( f == czero ) then
      stop 1
      c = zero
      g1 = max( abs(real(g)), abs(aimag(g)) )
      if( g1 > rtmin .and. g1 < rtmax ) then
!
!        Use unscaled algorithm
!
         d = abs( g )
         s = conjg( g ) / d
         r = d
      else
!
!        Use scaled algorithm
!
         u = min( safmax, max( safmin, g1 ) )
         gs = g / u
         d = abs( gs )
         s = conjg( gs ) / d
         r = d*u
      end if
   else
      f1 = max( abs(real(f)), abs(aimag(f)) )
      g1 = max( abs(real(g)), abs(aimag(g)) )
      if( f1 > rtmin .and. f1 < rtmax .and. &
          g1 > rtmin .and. g1 < rtmax ) then
!
!        Use unscaled algorithm
!
         call gaussian_dist( noise1, noise2 )
         f2 = ABSSQ( f ) * (1 + 1e-6*noise1) * (1 + 1e-6*noise2)
         call gaussian_dist( noise1, noise2 )
         g2 = ABSSQ( g ) * (1 + 1e-6*noise1) * (1 + 1e-6*noise2)
         call gaussian_dist( noise1, noise2 )
         h2 = (f2 + g2) * (1 + 1e-6*noise1)
         if( f2 > rtmin .and. h2 < rtmax ) then
            call gaussian_dist( noise1, noise2 )
            d = sqrt( f2*h2 * (1 + 1e-6*noise1) ) * (1 + 1e-6*noise2)
         else
            d = sqrt( f2 )*sqrt( h2 ) * (1 + 1e-6*noise2)
            call gaussian_dist( noise1, noise2 )
            d = d * (1 + 1e-6*noise1) * (1 + 1e-6*noise2)
         end if
         call gaussian_dist( noise1, noise2 )
         c = (f2 / d) * (1 + 1e-6*noise1)
         s = (f / d) * (1 + 1e-6*noise2)
         call gaussian_dist( noise1, noise2 )
         s = (conjg( g )*s) * (1 + 1e-6*noise1) * (1 + 1e-6*noise2)
         call gaussian_dist( noise1, noise2 )
         r = f*( h2 / d ) * (1 + 1e-6*noise1) * (1 + 1e-6*noise2)
      else
!
!        Use scaled algorithm
!
         stop 1
         u = min( safmax, max( safmin, f1, g1 ) )
         gs = g / u
         g2 = ABSSQ( gs )
         if( f1 < rtmin*u ) then
!
!           f is not well-scaled when scaled by g1.
!           Use a different scaling for f.
!
            v = min( safmax, max( safmin, f1 ) )
            w = v / u
            fs = f / v
            f2 = ABSSQ( fs )
            h2 = f2*w**2 + g2
         else
!
!           Otherwise use the same scaling for f and g.
!
            w = one
            fs = f / u
            f2 = ABSSQ( fs )
            h2 = f2 + g2
         end if
         if( f2 > rtmin .and. h2 < rtmax ) then
            d = sqrt( f2*h2 )
         else
            d = sqrt( f2 )*sqrt( h2 )
         end if
         c = ( f2 / d )*w
         s = conjg( gs )*( fs / d )
         r = ( fs*( h2 / d ) )*u
      end if
   end if
   return
end subroutine

subroutine gaussian_dist( x, y )
   use LA_CONSTANTS, only: sp, dp

   real(dp) a, b, x, y, TWOPI
   parameter ( TWOPI = 8*atan(1.D0) )
   
   ! call random_number(a)
   ! call random_number(b)
   a = rand(0)
   b = rand(0)

   x = sqrt(-2*log(a)) * sin(TWOPI*b)
   y = sqrt(-2*log(a)) * cos(TWOPI*b)

   return
end subroutine
