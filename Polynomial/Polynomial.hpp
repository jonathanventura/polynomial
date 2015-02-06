#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include <Eigen/Core>
#include <unsupported/Eigen/Polynomials>
#include <vector>

#include <Polynomial/PolynomialInternal.hpp>

namespace Polynomial
{
    
    template<int deg>
    class Polynomial
    {
    public:
        double coef[deg+1];
        
        Eigen::Map< Eigen::Matrix<double,1,deg+1> > coefvec()
        {
            return Internal::vecmap<deg+1>(coef);
        }
        
        Eigen::Map< const Eigen::Matrix<double,1,deg+1> > coefvec() const
        {
            return Internal::vecmap<deg+1>(coef);
        }
        
        template <int start,int length>
        Eigen::Map< Eigen::Matrix<double,1,length> > coefblock()
        {
            return Eigen::Map< Eigen::Matrix<double,1,length> >(coef+start);
        }
        
        template <int start,int length>
        const
        Eigen::Map< const Eigen::Matrix<double,1,length> > coefblock() const
        {
            return Eigen::Map< const Eigen::Matrix<double,1,length> >(coef+start);
        }
        
        Polynomial()
        {
            coefvec() = Eigen::Matrix<double,1,deg+1>::Zero();
        }
        
        Polynomial(const Eigen::Matrix<double,1,deg+1> &coefin)
        {
            coefvec() = coefin;
        }
        
        Polynomial(const Polynomial<deg> &polyin)
        {
            coefvec() = polyin.coefvec();
        }
        
        Polynomial(const double *coefin)
        {
            std::copy(coefin,coefin+deg+1,coef);
        }
        
        template<int degin>
        Polynomial<Internal::max<degin,deg>::value> operator+(const Polynomial<degin> &poly) const
        {
            Polynomial<Internal::max<degin,deg>::value> p;
            p.coefvec().tail(degin+1) = poly.coefvec();
            p.coefvec().tail(deg+1) += coefvec();
            return p;
        }
        
        template<int degin>
        Polynomial<Internal::max<degin,deg>::value> operator-(const Polynomial<degin> &poly) const
        {
            Polynomial<Internal::max<degin,deg>::value> p;
            p.coefvec().tail(deg+1) = coefvec();
            p.coefvec().tail(degin+1) -= poly.coefvec();
            return p;
        }
        
        template<int degin>
        Polynomial<degin+deg> operator*(const Polynomial<degin> &poly) const
        {
            Polynomial<degin+deg> p;
            Internal::PartialConv<deg,degin>::compute(p.coef,coef,poly.coef);
            return p;
        }
        
        Polynomial<deg> operator*(const double c) const
        {
            return Polynomial<deg>(coefvec()*c);
        }
        
        double eval(double x) const
        {
            return Internal::PolyVal<deg>::compute(coef,x);
        }
        
        void realRoots(std::vector<double> &roots) const
        {
            if ( coef[0] == 0 )
            {
                Eigen::PolynomialSolver<double,deg-1> ps;
                ps.compute(coefvec().tail(deg).reverse());
                ps.realRoots(roots);
            } else {
                Eigen::PolynomialSolver<double,deg> ps;
                ps.compute(coefvec().reverse());
                ps.realRoots(roots);
            }
        }
        
        void realRootsSturm(const double lb, const double ub, std::vector<double> &roots) const
        {
            if ( coef[0] == 0 )
            {
                Internal::SturmRootFinder<deg-1> sturm( coef+1 );
                sturm.realRoots( lb, ub, roots );
            } else {
                Internal::SturmRootFinder<deg> sturm( coef );
                sturm.realRoots( lb, ub, roots );
            }
        }
        
        void realRootBounds( double &lb, double &ub )
        {
        }
    };
    
} // end namespace Polynomial

#endif

