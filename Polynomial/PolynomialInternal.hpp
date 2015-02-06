#ifndef POLYNOMIAL_INTERNAL_HPP
#define POLYNOMIAL_INTERNAL_HPP

#include <Eigen/Core>
#include <vector>
#include <queue>

namespace Polynomial
{
    namespace Internal
    {
        
        template<int a,int b> struct max{ static const int value=(a<b)?b:a; };
        template<int a,int b> struct min{ static const int value=(a>b)?b:a; };
        
        template <int length>
        inline
        Eigen::Map< Eigen::Matrix<double,1,length> > vecmap(double *data)
        {
            return Eigen::Map< Eigen::Matrix<double,1,length> >(data);
        }
        
        template <int length>
        inline
        Eigen::Map< const Eigen::Matrix<double,1,length> > vecmap(const double *data)
        {
            return Eigen::Map< const Eigen::Matrix<double,1,length> >(data);
        }
        
        template <int deg1,int deg2,int k>
        struct
        _PartialConv
        {
            static void
            compute( double *coefout, const double *coef1, const double *coef2 )
            {
                coefout[k] =
                (
                 vecmap<min<k,deg1>::value-max<0,k-deg2>::value+1>
                 (
                  coef1+max<0,k-deg2>::value
                  )
                 .array()
                 *
                 vecmap<min<k,deg1>::value-max<0,k-deg2>::value+1>
                 (
                  coef2+k-min<k,deg1>::value
                  ).reverse().array()
                 ).sum();
                _PartialConv<deg1,deg2,k-1>::compute( coefout, coef1, coef2 );
            }
        };
        
        template <int deg1,int deg2>
        struct
        _PartialConv<deg1,deg2,-1>
        {
            static void
            compute( double *coefout, const double *coef1, const double *coef2 )
            {
            }
        };
        
        template <int deg1,int deg2>
        struct
        PartialConv
        {
            static void
            compute( double *coefout, const double *coef1, const double *coef2 )
            {
                _PartialConv<deg1,deg2,deg1+deg2>::compute(coefout,coef1,coef2);
            }
        };
        
        template <int deg,int k>
        struct
        _PolyVal
        {
            static double
            compute( const double *coef, const double &x )
            {
                return coef[k] + x*_PolyVal<deg,k-1>::compute(coef,x);
            }
        };
        
        template <int deg>
        struct
        _PolyVal<deg,0>
        {
            static double
            compute( const double *coef, const double &x )
            {
                return coef[0];
            }
        };
        
        template <int deg>
        struct
        PolyVal
        {
            static double
            compute( const double *coef, const double &x )
            {
                return Internal::_PolyVal<deg,deg>::compute( coef, x );
            }
        };
        
        template<int deg>
        struct BisectionMethod
        {
            static double compute(const double *coef, const double lower, const double upper )
            {
                double a = lower;
                double fa = PolyVal<deg>::compute(coef,a);
                char fapos = (fa>0);
                double b = upper;
                double c = 0.5*(a+b);
                
                for ( int i = 0; i < 24; i++ )
                {
                    const double fc = PolyVal<deg>::compute(coef,c);
                    const char fcpos = (fc>0);
                    if ( fcpos == fapos )
                    {
                        a = c;
                        fa = fc;
                    } else {
                        b = c;
                    }
                    c = 0.5*(a+b);
                }
                
                return c;
            }
        };
        
        template<int deg>
        struct SturmChain
        {
            double coef[deg+1];
            SturmChain<deg-1> next;
        };
        
        template<>
        struct SturmChain<0>
        {
            double coef[1];
        };
        
        template<int degin,int degout,int k>
        struct _GetSturmChain
        {
            static SturmChain<degout>& compute( SturmChain<degin> &sc )
            {
                return _GetSturmChain<degin-1,degout,k-1>::compute(sc.next);
            }
        };
        
        template<int degin,int degout>
        struct _GetSturmChain<degin,degout,0>
        {
            static SturmChain<degout>& compute( SturmChain<degin> &sc )
            {
                return sc;
            }
        };
        
        template<int degin,int degout>
        struct GetSturmChain
        {
            static SturmChain<degout>& compute( SturmChain<degin> &sc )
            {
                return _GetSturmChain<degin,degout,degin-degout>::compute(sc);
            }
        };
        
        template <int deg,int k>
        struct _ModPoly
        {
            static const void compute( SturmChain<deg> &sc )
            {
                const Eigen::Map< Eigen::Matrix<double,k+3,1> > u( GetSturmChain<deg,k+2>::compute(sc).coef );
                const Eigen::Map< Eigen::Matrix<double,k+2,1> > v( GetSturmChain<deg,k+1>::compute(sc).coef );
                Eigen::Matrix<double,k+3,1> r = u;
                const double s = v(0)/fabs(v(0));
                r(k+2) *= s;
                const int blocklength = k+1;
                r.block(1,0,blocklength,1) = s * r.block(1,0,blocklength,1) - r(0) * v.block(1,0,blocklength,1);
                r.block(2,0,blocklength,1) = s * r.block(2,0,blocklength,1) - r(1) * v.block(1,0,blocklength,1);
                Eigen::Map< Eigen::Matrix<double,k+1,1> > rout( GetSturmChain<deg,k>::compute(sc).coef );
                rout = r.block(2,0,k+1,1);
                double f = -fabs(rout(0));
                rout /= f;
                _ModPoly<deg,k-1>::compute(sc);
            }
        };
        
        template <int deg>
        struct _ModPoly<deg,0>
        {
            static const void compute( SturmChain<deg> &sc )
            {
                const Eigen::Map< Eigen::Matrix<double,3,1> > u( GetSturmChain<deg,2>::compute(sc).coef );
                const Eigen::Map< Eigen::Matrix<double,2,1> > v( GetSturmChain<deg,1>::compute(sc).coef );
                Eigen::Matrix<double,3,1> r = u;
                const double s = v(0)/fabs(v(0));
                r(2) *= s;
                r(1) = s * r(1) - r(0) * v(1);
                r(2) = s * r(2) - r(1) * v(1);
                GetSturmChain<deg,0>::compute(sc).coef[0] = -r(2);
            }
        };
        
        template <int deg>
        struct ModPoly
        {
            static const void compute( SturmChain<deg> &sc )
            {
                _ModPoly<deg,deg-2>::compute(sc);
            }
        };
        
        template <int deg>
        struct _SignChanges
        {
            static const void compute( SturmChain<deg> &sc, const double x, const double lf, int &count )
            {
                double f = PolyVal<deg>::compute( sc.coef, x );
                count += (f*lf<0);
                _SignChanges<deg-1>::compute( sc.next, x, f, count );
            }
        };
        
        template <>
        struct _SignChanges<0>
        {
            static const void compute( SturmChain<0> &sc, const double x, const double lf, int &count )
            {
                double f = sc.coef[0];
                count += (f*lf<0);
            }
        };
        
        template <int deg>
        struct SignChanges
        {
            static const void compute( SturmChain<deg> &sc, const double x, int &count )
            {
                double f = PolyVal<deg>::compute( sc.coef, x );
                count = 0;
                _SignChanges<deg-1>::compute( sc.next, x, f, count );
            }
        };
        
        template <>
        struct SignChanges<0>
        {
            static const void compute( SturmChain<0> &sc, const double x, int &count )
            {
                count = 0;
            }
        };
        
        template<int deg>
        class SturmRootFinder
        {
            Eigen::VectorXd modp( const Eigen::VectorXd &u, const Eigen::VectorXd &v )
            {
                int uord = u.rows()-1;
                Eigen::VectorXd r = u;
                double s = v(0)/fabs(v(0));
                r(uord) *= s;
                int blocklength = uord-1;
                r.block(1,0,blocklength,1) = s * r.block(1,0,blocklength,1) - r(0) * v.block(1,0,blocklength,1);
                r.block(2,0,blocklength,1) = s * r.block(2,0,blocklength,1) - r(1) * v.block(1,0,blocklength,1);
                return r.block(2,0,blocklength,1);
            }
            
            struct SturmInterval
            {
                double lb, ub;
                int sclb, scub;
                SturmInterval(const double _lb, const double _ub, const int _sclb, const int _scub)
                : lb(_lb), ub(_ub), sclb(_sclb), scub(_scub) { }
            };
            
            struct RootInterval
            {
                double lb, ub;
                RootInterval(const double _lb, const double _ub)
                : lb(_lb), ub(_ub) { }
            };
            
            int signchanges(double x)
            {
                double f[deg+1];
                for ( int i = 0; i < deg+1; i++ ) f[i] = PolyValDynamic(s[i],x);
                int n = 0;
                for ( int i = 1; i < deg+1; i++ ) n += (f[i]*f[i-1]<0);
                return n;
            }
            
            void bisection(const double lb, const double ub, std::vector<double> &roots)
            {
                int sclb = 0;
                int scub = 0;
                SignChanges<deg>::compute(sc,lb,sclb);
                SignChanges<deg>::compute(sc,ub,scub);
                std::queue<SturmInterval> intervals;
                intervals.push( SturmInterval(lb,ub,sclb,scub) );
                
                std::vector<RootInterval> root_intervals;
                
                while ( !intervals.empty() )
                {
                    SturmInterval interval = intervals.front();
                    intervals.pop();
                    
                    // get num roots in interval
                    int n = interval.sclb-interval.scub;
                    
                    if ( n == 0 )
                    {
                        
                    }
                    else if ( n == 1 )
                    {
                        root_intervals.push_back( RootInterval(interval.lb,interval.ub) );
                    }
                    else if ( n == 2 )
                    {
                        double m = (interval.lb+interval.ub)*0.5;
                        
                        // get sign changes at m
                        int scm = 0;
                        SignChanges<deg>::compute(sc,m,scm);
                        
                        // add [lb,m] interval
                        intervals.push( SturmInterval(interval.lb,m,interval.sclb,scm) );
                        
                        // add [m,ub] interval
                        intervals.push( SturmInterval(m,interval.ub,scm,interval.scub) );
                    }
                }
                
                roots.resize(root_intervals.size());
                for ( int i = 0; i < root_intervals.size(); i++ )
                {
                    roots[i] = BisectionMethod<deg>::compute(sc.coef,root_intervals[i].lb,root_intervals[i].ub);
                }
            }
        public:
            Eigen::VectorXd s[deg+1];
            SturmChain<deg> sc;
            
            SturmRootFinder( const double *coef )
            {
                std::copy( coef, coef+deg+1, sc.coef );
                std::copy( coef, coef+deg, sc.next.coef );
                double f = fabs(sc.coef[0] * deg);
                for ( int i = 0; i < deg; i++ ) sc.next.coef[i] *= (deg-i)/f;
                ModPoly<deg>::compute(sc);
            }
            
            void realRoots(const double lb, const double ub, std::vector<double> &roots)
            {
                bisection(lb,ub,roots);
            }
        };

    } // end namespace Internal
    
} // end namespace Polynomial

#endif
