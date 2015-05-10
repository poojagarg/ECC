#include "../src/ECC.h"

using namespace Givaro;
using namespace std;

int main()
{
	SignatureECDSA S;
	Key K;
	ellipticCurveFq E_Fq;//define Elliptic Curve over Fq
	ExtensionField::Element x,y;
	std::cout<<"give (x,y) of Base Point: ";
    E_Fq.field->readElement(x);
    E_Fq.field->readElement(y);
    ecPoint basePoint(x,y);
	Integer n,message;
	cin>>n;
	cin>>message;
	
	KeyPairGeneration(K, basePoint,E_Fq,n);//ord(P)=n, generates (Q,d such that Q=d*P)
	SignatureGeneration(S, K,message,basePoint,E_Fq,n);
	showSignature(S);
	inputSignature(S);

	if(SignatureVerification(S,message,basePoint,K.Q,E_Fq,n))
	{
		cout<<"Signature is valid";
	}
	else
	{
		cout<<"Signature is invalid";
	}
	
	return 0;
}