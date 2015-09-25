Polynomial
==========

Polynomial is a heavily templated C++ class for polynomials.  It provides functionality for polynomial algebra (addition, subtraction, and multiplication) as well as root-finding using either the numerically-accurate companion matrix method or the much faster Sturm sequences method.  It supports both statically-sized and dynamically-sized polynomials.

Wherever possible, Polynomial avoids for loops and instead uses template patterns for statically sized vectors to maximize speed.

Dependencies
------------

Polynomial depends on the [Eigen](http://eigen.tuxfamily.org/) linear algebra package.

Usage Example
=============

    #include <Polynomial/Polynomial.hpp>

    Eigen::Matrix<double,6,1> acoeffs;
    acoeffs << 1,2,3,4,5,6;
    Polynomial<5> a( acoeffs );
    
    Eigen::Matrix<double,5,1> bcoeffs;
    bcoeffs << 1,2,3,4,5;
    Polynomial<4> b( bcoeffs );

    Polynomial<9> c = a+b;
    std::vector<double> roots = c.realRoots();
    

