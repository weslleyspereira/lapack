/*
 * Copyright (c) 2020 Christoph Conrads (https://christoph-conrads.name)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of the copyright holders nor the
 *    names of its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef LAPACK_TESTS_COMPLEX_HPP
#define LAPACK_TESTS_COMPLEX_HPP

/**
 * This file replaces the C++ standard library header `complex` because prior
 * to GCC 8, this header included the C standard library header `complex.h`
 * which causes various issues because of the predefined macros, e.g., the
 * macro `I`. This in turn breaks compilation. Here is an example with GCC 7.5
 * on Ubuntu 18 (Bionic):
 *
 * ```
 * In file included from /usr/include/c++/7/complex.h:36:0,
 *                 from /home/dev/lapack/LAPACKE/include/lapack.h:55,
 *                 from /home/dev/lapack/LAPACKE/include/lapacke.h:37,
 *                 from /home/dev/lapack/TESTING/cpp/lapack.hpp:34,
 *                 from /home/dev/lapack/TESTING/cpp/xGGQRCS.hpp:31,
 *                 from /home/dev/lapack/TESTING/cpp/xGGQRCS_tests.cpp:30:
 * /usr/include/boost/operators.hpp:307:26: error: expected identifier before ‘(’ token
 * template <class T, class I, class R, class B = operators_detail::empty_base<T> >
 * ```
 * The C++ code will be expanded by the preprocessor to
 * ```
 * template <class T, class (__extension__ 1.0iF), class R, class B = operators_detail::empty_base<T> >
 * ```
 *
 * The sole purpose of this header is to prevent this problem.
 *
 * References:
 * * GCC PR 82417 "Macros from C99 <complex.h> defined in C++11" (https://gcc.gnu.org/PR82417).
 */

#include <complex>

#if (defined __GNUC__ && __GNUC__ <= 7)
// for some reason the fix does not work without explicitly including complex.h
#include <complex.h>
#undef complex
#undef I
#undef R
#endif

#endif
