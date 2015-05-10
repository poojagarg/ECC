#include "ECC.h"

using namespace std;

//Algorithm 4.24 Key pair generation
void convertElementToInt(Integer& result,ExtensionField::Element& R, Integer& mod, ellipticCurveFq& E_Fq)
{
	result=0;
	ExtensionField::Element P;
	if (R.size()) 
	{
        E_Fq.field->Fp_X.assign(P, R);
        E_Fq.field->Fp_X.setdegree(P);
        for(unsigned long l=0;l<P.size();++l) 
        {
    		result+=P[l]*(l+1);//if P[l] loses precision, that's not a problem
            result=result%mod;
        }                    
    }
}
Integer H(Integer& message,Integer& n)
{
	return message%n;
}

void KeyPairGeneration(Key& K, ecPoint& P,ellipticCurveFq& E_Fq,Integer& n)//ord(P)=n
{
	Integer::random_lessthan (K.d,n);
	E_Fq.scalarMultiply(K.Q,P,K.d,n);
}

/*Algorithm 4.25 Public Key validation
INPUT: Domain parameters D = (q,FR, S,a,b, P,n,h), public Key Q. 
OUTPUT: Acceptance or rejection of the validity of Q.*/

bool publicKeyValidation(ecPoint& Q,ellipticCurveFq& E_Fq,Integer& n)//Q belongs to <P>, ord(P)=n,
{
	if(Q.identity)
		return false;
	if(!(E_Fq.field->isElement(Q.x)&&E_Fq.field->isElement(Q.y)))
		return false;
	if(!E_Fq.verifyPoint(Q))
		return false;
	ecPoint T;
	E_Fq.scalarMultiply(T,Q,n,n);
	if(!T.identity)
		return false;
	return true;
}

/*Algorithm 4.29 ECDSA signature generation*/
//Assumption: H(x) gives a number of bitlength <= n
void SignatureGeneration(SignatureECDSA& S, Key& K,Integer& message,ecPoint& P,ellipticCurveFq& E_Fq,Integer& n)
{
	Integer k,_x=1,e;
	ecPoint temp;
	do
	{
		do
		{
			Integer::random_lessthan (k,n);	
			E_Fq.scalarMultiply(temp,P,k,n);
			convertElementToInt(_x,temp.x,n,E_Fq);
			S.r=_x%n;
		}while(S.r==0);
		e=H(message,n);
		invin(k,n);
		S.s=((e+K.d*S.r)*k)%n;
	}while(S.s==0);
}
/*Algorithm 4.30 ECDSA signature verification*/
bool SignatureVerification(SignatureECDSA& S,Integer& message,ecPoint& P,ecPoint& Q,ellipticCurveFq& E_Fq,Integer& n)
{
	if(S.r>n-1||S.r<1||S.s>n-1||S.s<1)
		return false;
	Integer e,w,u1,u2;
	ecPoint X;
	e=H(message,n);
	inv(w,S.s,n);
	u1=(e*w)%n;
	u2=(S.r*w)%n;
	E_Fq.compute(X,u1,P,u2,Q,n);
	if(X.identity)
		return false;
	convertElementToInt(e,X.x,n,E_Fq);
	w=e%n;	
	if(w==S.r)
		return true;
	return false;		
}
//Assumption: m<n
void SignatureGeneration(SignatureELGAMAL& S, Key& K,Integer& message,ecPoint& P,ellipticCurveFq& E_Fq,Integer& n)
{
	Integer k,_x;
	ecPoint temp;
	do
	{
		Integer::random_lessthan (k,n);	
	}while(gcd(k,n)!=1);
	E_Fq.scalarMultiply(S.R,P,k,n);
	convertElementToInt(_x,S.R.x,n,E_Fq);//x co-ordinate of point R of Signature S
	invin(k,n);
	S.s=((message-K.d*_x)*k)%n;
	if(S.s<=0)
		S.s+=n;
}
bool SignatureVerification(SignatureELGAMAL& S,Integer& message,ecPoint& P,ecPoint& Q,ellipticCurveFq& E_Fq,Integer& n)
{
	ecPoint V1,V2;
	Integer _x;
	convertElementToInt(_x,S.R.x,n,E_Fq);//x co-ordinate of point R of Signature S
	E_Fq.compute(V1,_x,Q,S.s,S.R,n);
	E_Fq.scalarMultiply(V2,P,message,n);
	if(V1==V2)
		return true;
	else
		return false;
}
void EncryptionELGAMAL(ecPoint& M1, ecPoint& M2, ecPoint& P,ecPoint& Q,ellipticCurveFq& E_Fq,Integer& n, ecPoint& message)
{
	Integer k;
	ecPoint temp;
	Integer::random_lessthan (k,n);
	E_Fq.scalarMultiply(M1,P,k,n);
	E_Fq.scalarMultiply(temp,Q,k,n);
	E_Fq.add(M2,message,temp);
}
void DecryptionELGAMAL(ecPoint& message, Integer& d, ecPoint& P,ecPoint& Q,ellipticCurveFq& E_Fq,Integer& n, ecPoint& M1, ecPoint& M2)
{
	ecPoint temp;
	E_Fq.scalarMultiply(temp,M1,d,n);
	E_Fq.inv(M1,temp);
	E_Fq.add(message,M2,M1);
}
/*Output: Ka- public Key with Alice, Kb- public Key with Bob
a and b: secret integer by Alice and Bob respectively*/
void DiffieHellmanKeyExchange(Integer& a, Integer& b, ecPoint& Ka,ecPoint& Kb,ecPoint& P,ellipticCurveFq& E_Fq,Integer& n)
{	
	
	ecPoint Pa,Pb;
	
	Integer::random_lessthan (a,n);
	E_Fq.scalarMultiply(Pa,P,a,n);

	Integer::random_lessthan (b,n);
	E_Fq.scalarMultiply(Pb,P,b,n);

	E_Fq.scalarMultiply(Ka,Pb,a,n);
	E_Fq.scalarMultiply(Kb,Pa,b,n);

}
void inputSignature(SignatureELGAMAL& S,ellipticCurveFq& E_Fq)
{
	cin>>S.s;
	std::cout<<"give (x,y) of Point R in Signature: ";
    E_Fq.field->readElement(S.R.x);
    E_Fq.field->readElement(S.R.y);
}
void inputSignature(SignatureECDSA& S)
{
	cin>>S.s>>S.r;
}
void showSignature(SignatureECDSA& S)
{
	cout<<"s: "<<S.s<<"\t"<<"r: "<<S.r;
}
void showSignature(SignatureELGAMAL& S,ellipticCurveFq& E_Fq)
{
	cout<<"s: "<<S.s<<endl;
	E_Fq.show(S.R);
}
/*assumption: 0 ≤ message < (p*m)/100, where E(Fp^m) */
/*void integerToPoint(ecPoint& M,Integer& message, ellipticCurveFq& E_Fq)
{
		Integer size;
		ExtensionField::Element Xj,acc,messageE,temp;
		pow(Integer& size, E_Fq->field.p, E_Fq->field.m);
		E_Fq->field.Fp_X.assign(acc,Fp_X.zero);
    ();

    	integerToElement(messageE,message,E_Fq->field);
		for(int j=0;j<100;j++)
		{
			E_Fq->field.scalarMultiply(temp,messageE,(Integer)100);	
			E_Fq->field.add(Xj,temp,acc);
			compute();
		}

}*/