#include "../src/ECC.h"

using namespace Givaro;
using namespace std;

int main()
{

	ecPoint M1, M2; 
	Key K;
	ellipticCurveFq E_Fq;//define Elliptic Curve over Fq
	ExtensionField::Element x,y;
	std::cout<<"give (x,y) of Base Point: ";
    E_Fq.field->readElement(x);
    E_Fq.field->readElement(y);
    ecPoint basePoint(x,y);
    std::cout<<"give (x,y) of Message: ";
    E_Fq.field->readElement(x);
    E_Fq.field->readElement(y);
    ecPoint message(x,y);
	Integer n;
	cin>>n;
	KeyPairGeneration(K, basePoint,E_Fq,n);//ord(basePoint)=n, generates (Q,d such that Q=d*basePoint)
	EncryptionELGAMAL(M1,M2,basePoint,K.Q,E_Fq,n,message);
	DecryptionELGAMAL(message,K.d,basePoint,K.Q,E_Fq,n,M1,M2);
	E_Fq.show(message);
	return 0;
}