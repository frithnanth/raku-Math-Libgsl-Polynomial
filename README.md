[![Build Status](https://travis-ci.org/frithnanth/raku-Math-Libgsl-Polynomial.svg?branch=master)](https://travis-ci.org/frithnanth/raku-Math-Libgsl-Polynomial)

NAME
====

Math::Libgsl::Polynomial - An interface to libgsl, the Gnu Scientific Library - Polynomials.

SYNOPSIS
========

```perl6
use Math::Libgsl::Raw::Polynomial :ALL;

use Math::Libgsl::Polynomial :ALL;
```

DESCRIPTION
===========

Math::Libgsl provides an interface to the polynomial evaluation in libgsl, the GNU Scientific Library.

Math::Libgsl::Polynomial makes these tags available:

  * :eval

  * :divdiff

  * :quad

  * :cubic

  * :complexsolve

### sub poly-eval(Positional $c, Num(Cool) $x --> Num) is export(:eval)

Evaluates a polynomial with real coefficients for the real variable x.

    my @c = 1, 2, 3;
    say poly-eval(@c, 10); # prints 321

### sub poly-eval-derivs(Num(Cool) $x, Int $maxn, Positional $c --> List) is export(:eval)

This function evaluates a polynomial and its derivatives. The output array contains the values of dⁿP(x)/dxⁿ for the specified value of x starting with n = 0.

### sub poly-dd(Positional $xa, Positional $ya, Positional $x --> List) is export(:divdiff)

This function computes a divided-difference representation of the interpolating polynomial for the points (xa, ya) and evaluates the polynomial for each point x.

### sub poly-dd-taylor(Positional $xa, Positional $ya, Num(Cool) $x --> List) is export(:divdiff)

This function converts the divided-difference representation of a polynomial to a Taylor expansion and evaluates the Taylor coefficients about the point x.

### sub poly-dd-hermite(Positional $xa, Positional $ya, Positional $dya, Positional $x --> List) is export(:divdiff)

This function computes a divided-difference representation of the interpolating Hermite polynomial for the points (xa, ya) and evaluates the polynomial for each point x.

### sub poly-solve-quadratic(Num(Cool) $a, Num(Cool) $b, Num(Cool) $c --> List) is export(:quad)

This function finds the real roots of the quadratic equation ax² + bx + c = 0. It returns a list of values: the number of the real roots found and zero, one or two roots; if present the roots are sorted in ascending order.

### sub poly-complex-solve-quadratic(Num(Cool) $a, Num(Cool) $b, Num(Cool) $c --> List) is export(:quad)

This function finds the complex roots of the quadratic equation ax² + bx + c = 0. It returns a list of values: the number of the real roots found and zero, one or two roots. The root are returned as Raku Complex values.

### sub poly-solve-cubic(Num(Cool) $a, Num(Cool) $b, Num(Cool) $c --> List) is export(:cubic)

This function finds the real roots of the cubic equation x³ + ax² + bx + c = 0. It returns a list of values: the number of the real roots found and one or three roots; the roots are sorted in ascending order.

### sub poly-complex-solve-cubic(Num(Cool) $a, Num(Cool) $b, Num(Cool) $c --> List) is export(:cubic)

This function finds the complex roots of the cubic equation x³ + ax² + bx + c = 0. The number of complex roots is returned (always three); the roots are returned in ascending order, sorted first by their real components and then by their imaginary components. The root are returned as Raku Complex values.

### sub poly-complex-solve(*@a --> List) is export(:complexsolve)

This function computes the roots of the general polynomial a₀ + a₁x + a₂x² + … + aₙ₋₁xⁿ⁻¹ The root are returned as Raku Complex values.

C Library Documentation
=======================

For more details on libgsl see [https://www.gnu.org/software/gsl/](https://www.gnu.org/software/gsl/). The excellent C Library manual is available here [https://www.gnu.org/software/gsl/doc/html/index.html](https://www.gnu.org/software/gsl/doc/html/index.html), or here [https://www.gnu.org/software/gsl/doc/latex/gsl-ref.pdf](https://www.gnu.org/software/gsl/doc/latex/gsl-ref.pdf) in PDF format.

Prerequisites
=============

This module requires the libgsl library to be installed. Please follow the instructions below based on your platform:

Debian Linux
------------

    sudo apt install libgsl23 libgsl-dev libgslcblas0

That command will install libgslcblas0 as well, since it's used by the GSL.

Ubuntu 18.04 and Ubuntu 20.04
------------

libgsl23 and libgslcblas0 have a missing symbol on Ubuntu 18.04. I solved the issue installing the Debian Buster version of those three libraries:

  * [http://http.us.debian.org/debian/pool/main/g/gsl/libgslcblas0_2.5+dfsg-6_amd64.deb](http://http.us.debian.org/debian/pool/main/g/gsl/libgslcblas0_2.5+dfsg-6_amd64.deb)

  * [http://http.us.debian.org/debian/pool/main/g/gsl/libgsl23_2.5+dfsg-6_amd64.deb](http://http.us.debian.org/debian/pool/main/g/gsl/libgsl23_2.5+dfsg-6_amd64.deb)

  * [http://http.us.debian.org/debian/pool/main/g/gsl/libgsl-dev_2.5+dfsg-6_amd64.deb](http://http.us.debian.org/debian/pool/main/g/gsl/libgsl-dev_2.5+dfsg-6_amd64.deb)

Installation
============

To install it using zef (a module management tool):

    $ zef install Math::Libgsl::Polynomial

AUTHOR
======

Fernando Santagata <nando.santagata@gmail.com>

COPYRIGHT AND LICENSE
=====================

Copyright 2020 Fernando Santagata

This library is free software; you can redistribute it and/or modify it under the Artistic License 2.0.

