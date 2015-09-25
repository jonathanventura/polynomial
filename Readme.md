Polynomial
==========

Polynomial is a heavily templated C++ class for polynomials.  It provides functionality for polynomial algebra (addition, subtraction, and multiplication) as well as root-finding using either the numerically-accurate companion matrix method or the much faster Sturm sequences method.  It supports both statically-sized and dynamically-sized polynomials.

Wherever possible, Polynomial avoids for loops and instead uses template patterns for statically sized vectors to maximize speed.

Dependencies
------------

Polynomial depends on the [Eigen](http://eigen.tuxfamily.org/) linear algebra package.

Install
-------

Polynomial is a pure template library so no compilation is needed.  However, a CMake file is provided for installing the headers and compiling the test program.

Usage Example
-------------

    #include <Polynomial/Polynomial.hpp>

    Eigen::Matrix<double,6,1> acoeffs;
    acoeffs << 1,2,3,4,5,6;
    Polynomial<5> a( acoeffs );
    
    Eigen::Matrix<double,5,1> bcoeffs;
    bcoeffs << 1,2,3,4,5;
    Polynomial<4> b( bcoeffs );

    Polynomial<9> c = a*b;
    std::vector<double> roots = c.realRoots();

Todo
----

Allow for adjustable precision (number of bisections) in Sturm sequences root finding.

Improve the computation of lower and upper root bounds.

Add support for lazy evaluation as in the Eigen library.

License
-------

Copyright (c) 2015, Jonathan Ventura

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

