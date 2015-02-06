#include <iostream>
#include <vector>

#include <Polynomial/Polynomial.hpp>

using Polynomial::Polynomial;

int main( int argc, char **argv )
{
/*
    Polynom<5> a;
    a.coefvec() << 1,2,3,4,5,6;
    Polynom<4> b;
    b.coefvec() << 1,2,3,4,5;

    std::cout << (a+b).coefvec() << "\n";
    std::cout << (a-b).coefvec() << "\n";
    
    std::cout << (a*b).coefvec() << "\n";

    std::cout << a.eval(5.f) << "\n";
*/
    Polynomial<20> a;
    FILE *f = fopen("/Users/jventura/code/solverlib/data/polydata.dat","r");
    int deg;
    fread(&deg,sizeof(int),1,f);
    fread(&a.coef,sizeof(double),21,f);
    fclose(f);

    std::vector<double> roots;
    a.realRoots(roots);
    std::cout << "companion matrix:\n";
    for ( int i = 0; i < roots.size(); i++ ) std::cout << roots[i] << "\n";

    roots.clear();
    a.realRootsSturm(-15*M_PI/180.,15*M_PI/180.,roots);
    std::cout << "sturm:\n";
    for ( int i = 0; i < roots.size(); i++ ) std::cout << roots[i] << "\n";

    return 0;
}
