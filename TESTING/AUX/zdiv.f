*> \brief zdiv tests the robustness and precision of the double complex division
*> \author Weslley S. Pereira, University of Colorado Denver, U.S.
*
*> \verbatim
*>
*> Real values for test:
*> (1) x = b/2, where b is the Blue's min constant. (b/2)^2 underflows but the norm is b/2.
*> (2) x = UN, where UN is the underflow threshold. UN^2 underflows but the norm is UN.
*> (3) x = SM, where SM is the smallest subnormal number. SM^2 underflows but the norm is SM.
*>     Mind that not all platforms might implement subnormal numbers.
*> (4) x = 2*B, where B is the Blue's max constant. (2*B)^2 overflows but the norm is 2*B.
*> (5) x = OV, where OV is the overflow threshold. OV^2 overflows but the norm is OV.
*>
*> Tests:
*> (a) y = x + 0 * I, y/y = 1
*> (b) y = 0 + x * I, y/y = 1
*> (c) y = x + x * I, y/y = 1
*> (d) y1 = 0 + x * I, y2 = x + 0 * I, y1/y2 = I
*> (e) y1 = 0 + x * I, y2 = x + 0 * I, y2/y1 = -I
*> (f) y = x + x * I, y/conj(y) = I
*>
*> Special cases:
*>
*> (g) Inf inputs:
*>    y = ( Inf + 0   * I)
*>    y = ( 0   + Inf * I)
*>    y = (-Inf + 0   * I)
*>    y = ( 0   - Inf * I)
*>    y = ( Inf + Inf * I)
*> Tests:
*>    0 / y is either 0 or NaN.
*>    1 / y is either 0 or NaN.
*>    y / y is NaN.
*>
*> (h) NaN inputs:
*>    y = (NaN + 0   * I)
*>    y = (0   + NaN * I)
*>    y = (NaN + NaN * I)
*> Tests:
*>    0 / y is NaN.
*>    1 / y is NaN.
*>    y / y is NaN.
*>
*> \endverbatim
*
      program zdiv

      integer           N, i, nNaN, nInf
      parameter       ( N = 5, nNaN = 3, nInf = 5 )

      double precision  X( N ), threeFourth, fiveFourth, aInf, aNaN
      parameter       ( threeFourth = 3.0d0 / 4,
     $                  fiveFourth = 5.0d0 / 4 )

      double complex    Y, Y2, R, cInf( nInf ), cNaN( nNaN ), czero,
     $                  cone
      parameter       ( czero = DCMPLX( 0.0d0, 0.0d0 ),
     $                  cone  = DCMPLX( 1.0d0, 0.0d0 ) )
*
      intrinsic         DCONJG, DBLE, RADIX, CEILING, TINY, DIGITS,
     $                  MAXEXPONENT, MINEXPONENT, FLOOR, HUGE, DCMPLX
*
      X(1) = DBLE(RADIX(0.0d0))**CEILING(
     $            (MINEXPONENT(0.0d0) - 1) * 0.5d0 )
      X(2) = TINY(0.0d0)
      X(3) = TINY(0.0d0) * DBLE(RADIX(0.0d0))**( 1-DIGITS(0.0d0) )
      X(4) = DBLE(RADIX(0.0d0))**FLOOR(
     $            (MAXEXPONENT(0.0d0) - DIGITS(0.0d0) + 1) * 0.5d0 )
      X(5) = HUGE(0.0d0)
*
      aInf = X(5) * 2
      cInf(1) = DCMPLX( aInf, 0.0d0 )
      cInf(2) = DCMPLX(-aInf, 0.0d0 )
      cInf(3) = DCMPLX( 0.0d0, aInf )
      cInf(4) = DCMPLX( 0.0d0,-aInf )
      cInf(5) = DCMPLX( aInf,  aInf )
*
      aNaN = aInf / aInf
      cNaN(1) = DCMPLX( aNaN, 0.0d0 )
      cNaN(2) = DCMPLX( 0.0d0, aNaN )
      cNaN(3) = DCMPLX( aNaN,  aNaN )
*
*     Test (a) y = x + 0 * I, y/y = 1
      do 10 i = 1, N
          Y = DCMPLX( X(i), 0.0d0 )
          R = Y / Y
          if( R .ne. 1.0D0 ) then
              WRITE( *, FMT = 9999 ) Y, Y, R, 1.0D0
          endif
  10  continue
*
*     Test (b) y = 0 + x * I, y/y = 1
      do 20 i = 1, N
          Y = DCMPLX( 0.0d0, X(i) )
          R = Y / Y
          if( R .ne. 1.0D0 ) then
              WRITE( *, FMT = 9999 ) Y, Y, R, 1.0D0
          endif
  20  continue
*
*     Test (c) y = x + x * I, y/y = 1
      do 30 i = 1, N
          Y = DCMPLX( X(i), X(i) )
          R = Y / Y
          if( R .ne. 1.0D0 ) then
              WRITE( *, FMT = 9999 ) Y, Y, R, 1.0D0
          endif
  30  continue
*
*     Test (d) y1 = 0 + x * I, y2 = x + 0 * I, y1/y2 = I
      do 40 i = 1, N
          Y  = DCMPLX( 0.0d0, X(i) )
          Y2 = DCMPLX( X(i), 0.0d0 )
          R = Y / Y2
          if( R .ne. DCMPLX(0.0D0,1.0D0) ) then
              WRITE( *, FMT = 9999 ) Y, Y2, R, DCMPLX(0.0D0,1.0D0)
          endif
  40  continue
*
*     Test (e) y1 = 0 + x * I, y2 = x + 0 * I, y2/y1 = -I
      do 50 i = 1, N
          Y  = DCMPLX( 0.0d0, X(i) )
          Y2 = DCMPLX( X(i), 0.0d0 )
          R = Y2 / Y
          if( R .ne. DCMPLX(0.0D0,-1.0D0) ) then
              WRITE( *, FMT = 9999 ) Y2, Y, R, DCMPLX(0.0D0,-1.0D0)
          endif
  50  continue
*
*     Test (f) y = x + x * I, y/conj(y) = I
      do 60 i = 1, N
          Y  = DCMPLX( X(i), X(i) )
          R = Y / DCONJG( Y )
          if( R .ne. DCMPLX(0.0D0,1.0D0) ) then
              WRITE( *, FMT = 9999 ) Y, DCONJG( Y ), R,
     $                               DCMPLX(0.0D0,1.0D0)
          endif
  60  continue
*
*     Test (g) Infs
      do 70 i = 1, nInf
          Y = cInf(i)
          R = czero / Y
          if( (R .ne. czero) .and. (R .eq. R) ) then
              WRITE( *, FMT = 9998 ) czero, Y, R
          endif
          R = cone / Y
          if( (R .ne. 0.0d0) .and. (R .eq. R) ) then
              WRITE( *, FMT = 9998 ) cone, Y, R
          endif
          R = Y / Y
          if( R .eq. R ) then
              WRITE( *, FMT = 9998 ) Y, Y, R
          endif
  70  continue
*
*     Test (h) NaNs
      do 80 i = 1, nNaN
          Y = cNaN(i)
          R = czero / Y
          if( R .eq. R ) then
              WRITE( *, FMT = 9998 ) czero, Y, R
          endif
          R = cone / Y
          if( R .eq. R ) then
              WRITE( *, FMT = 9998 ) cone, Y, R
          endif
          R = Y / Y
          if( R .eq. R ) then
              WRITE( *, FMT = 9998 ) Y, Y, R
          endif
  80  continue
*
 9998 FORMAT( ' ** (', (ES10.3,SP,ES10.3,"*I"), ' ) / ( ',
     $                 (ES10.3,SP,ES10.3,"*I"), ' ) = ',
     $                 (ES10.3,SP,ES10.3,"*I"), ' is unexpected' )
*
 9999 FORMAT( ' ** (', (ES10.3,SP,ES10.3,"*I"), ' ) / ( ',
     $                 (ES10.3,SP,ES10.3,"*I"), ' ) = ',
     $                 (ES10.3,SP,ES10.3,"*I"), ' differs from ',
     $                 (ES10.3,SP,ES10.3,"*I") )
*
*     End of zdiv
*
      END