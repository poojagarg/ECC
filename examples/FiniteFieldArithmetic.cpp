#include "../src/ellipticCurve.h"

using namespace Givaro;

int main()
{
    int p,m;
    std::cin>>p>>m;
    ExtensionField E(p, m);
    typedef ExtensionField::Element e;//polynomial element
    e A,B,R;
    std::cout<<"1-Add, 2-Sub, 3-Mul, 4-Div, 5-Inv, 6-Neg, 7-Sqr, 8-Scalar multiply";
    
    std::cout<<"Write A";
    E.readElement(A,true);
    std::cout<<"Write B";
    E.readElement(B,true);
    Integer n;
    std::cin>>n;
    switch(1)
    {
    	case 1:
    		E.add(R,A,B);
		    E.writeElement(A);
            std::cout<<" + ";
            E.writeElement(B);
            std::cout<<"=";
            E.writeElement(R);
            std::cout<<std::endl;
    	case 2:
    		E.sub(R,A,B);
		    E.writeElement(A);
            std::cout<<" - ";
            E.writeElement(B);
            std::cout<<"=";
            E.writeElement(R);
            std::cout<<std::endl;
    	case 3:
    		E.mul(R,A,B);
		    E.writeElement(A);
            std::cout<<" * ";
            E.writeElement(B);
            std::cout<<"=";
            E.writeElement(R);
            std::cout<<std::endl;
    	case 4:
    		E.div(R,A,B);
		    E.writeElement(A);
            std::cout<<" / ";
            E.writeElement(B);
            std::cout<<"=";
            E.writeElement(R);
            std::cout<<std::endl;
    	case 5:
    		E.inv(R,A);
            std::cout<<"inverse of ";
		    E.writeElement(A);
            std::cout<<"=";
            E.writeElement(R);
            std::cout<<std::endl;
    	case 6:
    		E.neg(R,A);
            std::cout<<"negation of ";
		    E.writeElement(A);
            std::cout<<"=";
            E.writeElement(R);
            std::cout<<std::endl;
    	case 7:
    		E.sqr(R,A);
		    std::cout<<"square of ";
            E.writeElement(A);
            std::cout<<"=";
            E.writeElement(R);
            std::cout<<std::endl;
    	case 8:
    		E.scalarMultiply(R,A,n);
		    std::cout<<n<<"(*)";
            E.writeElement(A);
            std::cout<<"=";
            E.writeElement(R);
            std::cout<<std::endl;
    }
}