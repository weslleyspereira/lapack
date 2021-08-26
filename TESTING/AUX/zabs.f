*  (1) x = b/2, where b is the Blue's min constant. (b/2)^2 underflows but the norm is b/2.
*  (2) x = UN, where UN is the underflow threshold. UN^2 underflows but the norm is UN.
*  (3) x = SM, where SM is the smallest subnormal number. SM^2 underflows but the norm is SM.
*      Mind that not all platforms might implement subnormal numbers.
*  (4) x = 2*B, where B is the Blue's max constant. (2*B)^2 overflows but the norm is 2*B.
*  (5) x = OV, where OV is the overflow threshold. OV^2 overflows but the norm is OV.
*
*  (a) y = x + 0 * I, |y| = x
*  (b) y = 0 + x * I, |y| = x
*  (c) y = (3/4)*x + x * I, |y| = x for cases (1),(2), and (4)
*
      program zabs

      integer           N, i, idxTestC(3), idx
      parameter       ( N = 5, idxTestC = (/1,2,4/) )

      double precision  X( N ), R, threeFourth, fiveFourth
      parameter       ( threeFourth = 3.0d0 / 4,
     $                  fiveFourth = 5.0d0 / 4 )

      double complex    Y
      intrinsic         abs, dble, radix, ceiling, tiny, digits,
     $                  maxexponent, minexponent, floor, huge, DCMPLX

      X(1) = dble(radix(0.0d0))**ceiling(
     $            (minexponent(0.0d0) - 1) * 0.5d0 )
      X(2) = tiny(0.0d0)
      X(3) = tiny(0.0d0) * dble(radix(0.0d0))**( 1-digits(0.0d0) )
      X(4) = dble(radix(0.0d0))**floor(
     $            (maxexponent(0.0d0) - digits(0.0d0) + 1) * 0.5d0 )
      X(5) = huge(0.0d0)

      print *, X

      do 10 i = 1, N
          Y = DCMPLX( X(i), 0.0d0 )
          R = abs( Y )
          print *, "abs( ", Y, " ) = ", R, " # should be ", X(i)
  10  continue

      do 20 i = 1, N
          Y = DCMPLX( 0.0d0, X(i) )
          R = abs( Y )
          print *, "abs( ", Y, " ) = ", R, " # should be ", X(i)
  20  continue
 
      do 30 idx = 1, 3
          i = idxTestC( idx )
          R = threeFourth * X(i)
          Y = DCMPLX( threeFourth * X(i), X(i) )
          R = abs( Y )
          print *, "abs( ", Y, " ) = ", R, " # should be ",
     $            fiveFourth*X(i)
  30  continue
*
*     End of zabs
*
      END