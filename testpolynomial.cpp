#include <iostream>
#include <vector>

#include <Polynomial/Polynomial.hpp>

using Polynomial::Polynomial;

int main( int argc, char **argv )
{
    Eigen::Matrix<double,6,1> acoeffs;
    acoeffs << 1,2,3,4,5,6;
    Polynomial<5> a( acoeffs );
    Eigen::Matrix<double,5,1> bcoeffs;
    bcoeffs << 1,2,3,4,5;
    Polynomial<4> b( bcoeffs );

    std::cout << (a+b).coefficients().transpose() << "\n";
    std::cout << (a-b).coefficients().transpose() << "\n";
    
    std::cout << (a*b).coefficients().transpose() << "\n";

    std::cout << a.eval(5.f) << "\n";
    
    
    Eigen::VectorXd adcoeffs(6);
    adcoeffs << 1,2,3,4,5,6;
    Polynomial<Eigen::Dynamic> ad( adcoeffs );
    Eigen::VectorXd bdcoeffs(5);
    bdcoeffs << 1,2,3,4,5;
    Polynomial<Eigen::Dynamic> bd( bdcoeffs );
    
    std::cout << (ad+bd).coefficients().transpose() << "\n";
    std::cout << (ad-bd).coefficients().transpose() << "\n";
    
    std::cout << (ad*bd).coefficients().transpose() << "\n";
    
    std::cout << ad.eval(5.f) << "\n";
    
    Eigen::Matrix<double,5,1> ccoeffs;
    ccoeffs <<    -0.8049,    -0.4430,    0.0938,    0.9150,    0.9298;
    Polynomial<4> c(ccoeffs);
    
    std::vector<double> roots;
    c.realRoots(roots);
    std::cout << "companion matrix:\n";
    for ( int i = 0; i < roots.size(); i++ ) std::cout << roots[i] << "\n";

    double lb,ub;
    c.rootBounds(lb,ub);
    std::cout << "root bounds: " << lb << "  " << ub << "\n";

    roots.clear();
    c.realRootsSturm(lb,ub,roots);
    std::cout << "sturm:\n";
    for ( int i = 0; i < roots.size(); i++ ) std::cout << roots[i] << "\n";

    exit(0);
    
    Eigen::Matrix<double,21,1> coef;
    
    FILE *f = fopen("/Users/jventura/code/solverlib/data/polydata.dat","r");
    int deg;
    fread(&deg,sizeof(int),1,f);
    fread(coef.data(),sizeof(double),21,f);
    fclose(f);
    Polynomial<20> poly( coef );
    
    Eigen::VectorXd coefd( coef );
    Polynomial<Eigen::Dynamic> polyd( coef );

    poly.rootBounds(lb,ub);
    std::cout << "root bounds: " << lb << "  " << ub << "\n";
    
    roots.clear();
    poly.realRoots(roots);
    std::cout << "companion matrix:\n";
    for ( int i = 0; i < roots.size(); i++ ) std::cout << roots[i] << "\n";

    roots.clear();
    poly.realRootsSturm(lb,ub,roots);
    std::cout << "sturm:\n";
    for ( int i = 0; i < roots.size(); i++ ) std::cout << roots[i] << "\n";

    double lbd,ubd;
    polyd.rootBounds(lbd,ubd);
    std::cout << "root bounds (dynamic): " << lbd << "  " << ubd << "\n";
    
    roots.clear();
    polyd.realRoots(roots);
    std::cout << "companion matrix (dynamic):\n";
    for ( int i = 0; i < roots.size(); i++ ) std::cout << roots[i] << "\n";
    
    roots.clear();
    polyd.realRootsSturm(lbd,ubd,roots);
    std::cout << "sturm (dynamic):\n";
    for ( int i = 0; i < roots.size(); i++ ) std::cout << roots[i] << "\n";

    return 0;
}
