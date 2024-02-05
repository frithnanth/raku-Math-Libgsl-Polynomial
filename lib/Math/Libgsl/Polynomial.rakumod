use v6;

unit module Math::Libgsl::Polynomial:ver<0.0.4>:auth<zef:FRITH>;

use Math::Libgsl::Raw::Polynomial :ALL;
use Math::Libgsl::Exception;
use Math::Libgsl::Raw::Complex;
use NativeCall;

# Polynomial evaluation
sub poly-eval(Positional $c, Num(Cool) $x --> Num) is export(:eval) {
  my CArray[num64] $in .= new: $c».Num;
  gsl_poly_eval($in, $c.elems, $x)
}
sub poly-eval-derivs(Num(Cool) $x, Int $maxk, Positional $c --> List) is export(:eval) {
  my CArray[num64] $in  .= new: $c».Num;
  my $out = CArray[num64].allocate: $maxk;
  my $res = gsl_poly_eval_derivs($in, $in.elems, $x, $out, $maxk);
  fail X::Libgsl.new: errno => $res if $res > 0;
  return $out.list;
}
# Divided difference representation of polynomials
sub poly-dd(Positional $xa, Positional $ya, Positional $x --> List) is export(:divdiff) {
  my CArray[num64] $inx .= new: $xa».Num;
  my CArray[num64] $iny .= new: $ya».Num;
  my $dd = CArray[num64].allocate: $inx.elems;
  my $res = gsl_poly_dd_init($dd, $inx, $iny, $inx.elems);
  fail X::Libgsl.new: errno => $res if $res > 0;
  my @out;
  for $x.list».Num -> $val {
    @out.push: gsl_poly_dd_eval($dd, $inx, $inx.elems, $val);
  }
  return @out;
}
sub poly-dd-taylor(Positional $xa, Positional $ya, Num(Cool) $x --> List) is export(:divdiff) {
  my CArray[num64] $inx .= new: $xa».Num;
  my CArray[num64] $iny .= new: $ya».Num;
  my $dd = CArray[num64].allocate: $inx.elems;
  my $res = gsl_poly_dd_init($dd, $inx, $iny, $inx.elems);
  fail X::Libgsl.new: errno => $res if $res > 0;
  my $c = CArray[num64].allocate: $inx.elems;
  my $w = CArray[num64].allocate: $inx.elems;
  $res = gsl_poly_dd_taylor($c, $x, $dd, $inx, $inx.elems, $w);
  fail X::Libgsl.new: errno => $res if $res > 0;
  return $c.list;
}
sub poly-dd-hermite(Positional $xa, Positional $ya, Positional $dya, Positional $x --> List) is export(:divdiff) {
  my CArray[num64] $inx  .= new: $xa».Num;
  my CArray[num64] $iny  .= new: $ya».Num;
  my CArray[num64] $indy .= new: $dya».Num;
  my $dd = CArray[num64].allocate: $xa.elems * 2;
  my $za = CArray[num64].allocate: $xa.elems * 2;
  my $res = gsl_poly_dd_hermite_init($dd, $za, $inx, $iny, $indy, $inx.elems);
  fail X::Libgsl.new: errno => $res if $res > 0;
  my @out;
  for $x.list».Num -> $val {
    @out.push: gsl_poly_dd_eval($dd, $za, $za.elems, $val);
  }
  return @out;
}
# Quadratic equations
sub poly-solve-quadratic(Num(Cool) $a, Num(Cool) $b, Num(Cool) $c --> List) is export(:quad) {
  my num64 $x0;
  my num64 $x1;
  my $ret = gsl_poly_solve_quadratic($a, $b, $c, $x0, $x1);
  return $ret, $x0, $x1;
}
sub poly-complex-solve-quadratic(Num(Cool) $a, Num(Cool) $b, Num(Cool) $c --> List) is export(:quad) {
  my gsl_complex $z0 = alloc_gsl_complex;
  my gsl_complex $z1 = alloc_gsl_complex;
  my $ret = gsl_poly_complex_solve_quadratic($a, $b, $c, $z0, $z1);
  my Complex $root0 = $z0.dat[0] + i * $z0.dat[1];
  my Complex $root1 = $z1.dat[0] + i * $z1.dat[1];
  return $ret, $root0, $root1;
}
# Cubic equations
sub poly-solve-cubic(Num(Cool) $a, Num(Cool) $b, Num(Cool) $c --> List) is export(:cubic) {
  my num64 $x0;
  my num64 $x1;
  my num64 $x2;
  my $ret = gsl_poly_solve_cubic($a, $b, $c, $x0, $x1, $x2);
  return $ret, $x0, $x1, $x2;
}
sub poly-complex-solve-cubic(Num(Cool) $a, Num(Cool) $b, Num(Cool) $c --> List) is export(:cubic) {
  my gsl_complex $z0 = alloc_gsl_complex;
  my gsl_complex $z1 = alloc_gsl_complex;
  my gsl_complex $z2 = alloc_gsl_complex;
  my $ret = gsl_poly_complex_solve_cubic($a, $b, $c, $z0, $z1, $z2);
  my Complex $root0 = $z0.dat[0] + i * $z0.dat[1];
  my Complex $root1 = $z1.dat[0] + i * $z1.dat[1];
  my Complex $root2 = $z2.dat[0] + i * $z2.dat[1];
  return $ret, $root0, $root1, $root2;
}
# General polynomial equations
sub poly-complex-solve(*@a --> List) is export(:complexsolve) {
  my CArray[num64] $ina .= new: @a».Num;
  my $z = CArray[num64].allocate(($ina.elems - 1) * 2);
  my gsl_poly_complex_workspace $w = gsl_poly_complex_workspace_alloc($ina.elems);
  my $ret = gsl_poly_complex_solve($ina, $ina.elems, $w, $z);
  fail X::Libgsl.new: errno => $ret if $ret > 0;
  my @out;
  for $z.list -> $re, $im {
    @out.push: $re + i * $im;
  }
  return @out;
}

=begin pod

=head1 NAME

Math::Libgsl::Polynomial - An interface to libgsl, the Gnu Scientific Library - Polynomials.

=head1 SYNOPSIS

=begin code :lang<raku>

use Math::Libgsl::Raw::Polynomial :ALL;

use Math::Libgsl::Polynomial :ALL;

=end code

=head1 DESCRIPTION

Math::Libgsl provides an interface to polynomial evaluation in libgsl, the GNU Scientific Library.

Math::Libgsl::Polynomial makes these tags available:

=item :eval
=item :divdiff
=item :quad
=item :cubic
=item :complexsolve

=head3 sub poly-eval(Positional $c, Num(Cool) $x --> Num) is export(:eval)

Evaluates a polynomial with real coefficients for the real variable x.
=begin code
my @c = 1, 2, 3;
say poly-eval(@c, 10); # prints 321
=end code

=head3 sub poly-eval-derivs(Num(Cool) $x, Int $maxn, Positional $c --> List) is export(:eval)

This function evaluates a polynomial and its derivatives.
The output array contains the values of dⁿP(x)/dxⁿ for the specified value of x starting with n = 0.

=head3 sub poly-dd(Positional $xa, Positional $ya, Positional $x --> List) is export(:divdiff)

This function computes a divided-difference representation of the interpolating polynomial for the points (xa, ya) and evaluates the polynomial for each point x.

=head3 sub poly-dd-taylor(Positional $xa, Positional $ya, Num(Cool) $x --> List) is export(:divdiff)

This function converts the divided-difference representation of a polynomial to a Taylor expansion and evaluates the Taylor coefficients about the point x.

=head3 sub poly-dd-hermite(Positional $xa, Positional $ya, Positional $dya, Positional $x --> List) is export(:divdiff)

This function computes a divided-difference representation of the interpolating Hermite polynomial for the points (xa, ya) and evaluates the polynomial for each point x.

=head3 sub poly-solve-quadratic(Num(Cool) $a, Num(Cool) $b, Num(Cool) $c --> List) is export(:quad)

This function finds the real roots of the quadratic equation ax² + bx + c = 0.
It returns a list of values: the number of the real roots found and zero, one or two roots; if present the roots are sorted in ascending order.

=head3 sub poly-complex-solve-quadratic(Num(Cool) $a, Num(Cool) $b, Num(Cool) $c --> List) is export(:quad)

This function finds the complex roots of the quadratic equation ax² + bx + c = 0.
It returns a list of values: the number of the real roots found and zero, one or two roots.
The root are returned as Raku Complex values.

=head3 sub poly-solve-cubic(Num(Cool) $a, Num(Cool) $b, Num(Cool) $c --> List) is export(:cubic)

This function finds the real roots of the cubic equation x³ + ax² + bx + c = 0.
It returns a list of values: the number of the real roots found and one or three roots; the roots are sorted in ascending order.

=head3 sub poly-complex-solve-cubic(Num(Cool) $a, Num(Cool) $b, Num(Cool) $c --> List) is export(:cubic)

This function finds the complex roots of the cubic equation x³ + ax² + bx + c = 0.
The number of complex roots is returned (always three); the roots are returned in ascending order, sorted first by their real components and then by their imaginary components.
The root are returned as Raku Complex values.

=head3 sub poly-complex-solve(*@a --> List) is export(:complexsolve)

This function computes the roots of the general polynomial a₀ + a₁x + a₂x² + … + aₙ₋₁xⁿ⁻¹
The root are returned as Raku Complex values.

=head1 C Library Documentation

For more details on libgsl see L<https://www.gnu.org/software/gsl/>.
The excellent C Library manual is available here L<https://www.gnu.org/software/gsl/doc/html/index.html>, or here L<https://www.gnu.org/software/gsl/doc/latex/gsl-ref.pdf> in PDF format.

=head1 Prerequisites

This module requires the libgsl library to be installed. Please follow the instructions below based on your platform:

=head2 Debian Linux and Ubuntu 20.04

=begin code
sudo apt install libgsl23 libgsl-dev libgslcblas0
=end code

That command will install libgslcblas0 as well, since it's used by the GSL.

=head2 Ubuntu 18.04

libgsl23 and libgslcblas0 have a missing symbol on Ubuntu 18.04.
I solved the issue installing the Debian Buster version of those three libraries:

=item L<http://http.us.debian.org/debian/pool/main/g/gsl/libgslcblas0_2.5+dfsg-6_amd64.deb>
=item L<http://http.us.debian.org/debian/pool/main/g/gsl/libgsl23_2.5+dfsg-6_amd64.deb>
=item L<http://http.us.debian.org/debian/pool/main/g/gsl/libgsl-dev_2.5+dfsg-6_amd64.deb>

=head1 Installation

To install it using zef (a module management tool):

=begin code
$ zef install Math::Libgsl::Polynomial
=end code

=head1 AUTHOR

Fernando Santagata <nando.santagata@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright 2020 Fernando Santagata

This library is free software; you can redistribute it and/or modify it under the Artistic License 2.0.

=end pod
