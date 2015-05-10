#include "../src/attacksECC.h"

using namespace std;

int main()
{
    //takes input parameter for E(Fq), where E is elliptic curve defined over K and Fq is extension of K
    ellipticCurveFq E_Fq;
    E_Fq.show();
    Integer n,result;
    ExtensionField::Element x1,y1,x2,y2,x3,y3;
    std::cout<<"Write x1";
    E_Fq.field->readElement(x1);
    std::cout<<"Write y1";
    E_Fq.field->readElement(y1);
    std::cout<<"Write x2";
    E_Fq.field->readElement(x2);
    std::cout<<"Write y2";
    E_Fq.field->readElement(y2);
    ecPoint P(x1,y1);
    ecPoint Q(x2,y2);
    E_Fq.show(P);
    E_Fq.show(Q);
    std::cout<<"order of P?"; 
    std::cin>>n;
    pohligHellman(result, P, Q,n,E_Fq);
    //if n is prime, can also use Pollard Rho method- pollardRho(result, P, Q,n,E_Fq);
    
    return 0;
}