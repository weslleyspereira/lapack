*  (1) x = b/2, where b is the Blue's min constant. (b/2)^2 underflows but the norm is b/2.
*  (2) x = UN, where UN is the underflow threshold. UN^2 underflows but the norm is UN.
*  (3) x = SM, where SM is the smallest subnormal number. SM^2 underflows but the norm is SM.
*      Mind that not all platforms might implement subnormal numbers.
*  (4) x = 2*B, where B is the Blue's max constant. (2*B)^2 overflows but the norm is 2*B.
*  (5) x = OV, where OV is the overflow threshold. OV^2 overflows but the norm is OV.
*
*  (a) y = x + 0 * I, y/y = 1
*  (b) y = 0 + x * I, y/y = 1
*  (c) y = x + x * I, y/y = 1
*  (d) y1 = 0 + x * I, y2 = 0 + x * I, y1/y2 = I
*  (e) y1 = 0 + x * I, y2 = 0 + x * I, y2/y1 = -I
*  (f) y = x + x * I, y/conj(y) = I
*
      program zdiv

      integer           N, i
      parameter       ( N = 5 )

      double precision  X( N ), threeFourth, fiveFourth
      parameter       ( threeFourth = 3.0d0 / 4,
     $                  fiveFourth = 5.0d0 / 4 )

      double complex    Y, Y2, R
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
              print *, Y, R
          endif
  20  continue
*
*     Test (c) y = x + x * I, y/y = 1
      do 30 i = 1, N
          Y = DCMPLX( X(i), X(i) )
          R = Y / Y
          if( R .ne. 1.0D0 ) then
              WRITE( *, FMT = 9999 ) Y, Y, R, 1.0D0
              print *, Y, R
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
              print *, Y, Y2, R
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
              print *, Y, Y2, R
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
              print *, Y, R
          endif
  60  continue
*
 9999 FORMAT( ' ** (', (ES10.3,SP,ES10.3,"*I"), ' ) / ( ',
     $                 (ES10.3,SP,ES10.3,"*I"), ' ) = ',
     $                 (ES10.3,SP,ES10.3,"*I"), ' differs from ',
     $                 (ES10.3,SP,ES10.3,"*I") )
*
*     End of zdiv
*
      END