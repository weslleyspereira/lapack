*  (1) x = b/2, where b is the Blue's min constant. (b/2)^2 underflows but the norm is b/2.
*  (2) x = UN, where UN is the underflow threshold. UN^2 underflows but the norm is UN.
*  (3) x = SM, where SM is the smallest subnormal number. SM^2 underflows but the norm is SM.
*      Mind that not all platforms might implement subnormal numbers.
*  (4) x = 2*B, where B is the Blue's max constant. (2*B)^2 overflows but the norm is 2*B.
*  (5) x = OV, where OV is the overflow threshold. OV^2 overflows but the norm is OV.
*
*  (a) y = x + 0 * I, |y| = x
*  (b) y = 0 + x * I, |y| = x
*  (c) y = (3/4)*x + x * I, |y| = (5/4)*x whenever (3/4)*x and (5/4)*x can be exactly stored
*  (d) y = (1/2)*x + (1/2)*x * I, |y| = (1/2)*x*sqrt(2) whenever (1/2)*x can be exactly stored
*
      program zabs

      integer           N, i
      parameter       ( N = 5 )

      double precision  X( N ), R, threeFourth, fiveFourth, answerC(N),
     $                  answerD(N), oneHalf
      parameter       ( threeFourth = 3.0d0 / 4,
     $                  fiveFourth = 5.0d0 / 4,
     $                  oneHalf = 1.0d0 / 2 )

      double complex    Y
      intrinsic         ABS, DBLE, RADIX, CEILING, TINY, DIGITS, SQRT,
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
      answerC(1) = fiveFourth * X(1)
      answerC(2) = fiveFourth * X(2)
      answerC(3) = -1
      answerC(4) = fiveFourth * X(4)
      answerC(5) = -1
*
      answerD(1) = (oneHalf * X(1)) * SQRT(2.0d0)
      answerD(2) = (oneHalf * X(2)) * SQRT(2.0d0)
      answerD(3) = -1
      answerD(4) = (oneHalf * X(4)) * SQRT(2.0d0)
      answerD(5) = (oneHalf * X(5)) * SQRT(2.0d0)
*
*     Test (a) y = x + 0 * I, |y| = x
      do 10 i = 1, N
          Y = DCMPLX( X(i), 0.0d0 )
          R = ABS( Y )
          if( R .ne. X(i) ) then
              WRITE( *, FMT = 9999 ) Y, R, X(i)
          endif
  10  continue
*
*     Test (b) y = 0 + x * I, |y| = x
      do 20 i = 1, N
          Y = DCMPLX( 0.0d0, X(i) )
          R = ABS( Y )
          if( R .ne. X(i) ) then
              WRITE( *, FMT = 9999 ) Y, R, X(i)
          endif
  20  continue
*
*     Test (c) y = (3/4)*x + x * I, |y| = (5/4)*x
      do 30 i = 1, N
          if( answerC(i) .lt. 0 ) go to 30
          Y = DCMPLX( threeFourth * X(i), X(i) )
          R = ABS( Y )
          if( R .ne. answerC(i) ) then
              WRITE( *, FMT = 9999 ) Y, R, answerC(i)
          endif
  30  continue
*
*     Test (d) y = (1/2)*x + (1/2)*x * I, |y| = (1/2)*x*sqrt(2)
      do 40 i = 1, N
          if( answerD(i) .lt. 0 ) go to 40
          Y = DCMPLX( oneHalf * X(i), oneHalf * X(i) )
          R = ABS( Y )
          if( R .ne. answerD(i) ) then
              WRITE( *, FMT = 9999 ) Y, R, answerD(i)
          endif
  40  continue
*
 9999 FORMAT( ' ** ABS( ', (ES10.3,SP,ES10.3,"*I"), ' ) = ', ES10.3,
     $        ' differs from ', ES10.3 )
*
*     End of zabs
*
      END