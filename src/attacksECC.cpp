/*INPUT: P ∈ E(Fq) of prime order n, Q ∈ ⟨P⟩. 
OUTPUT: The discrete logarithm l = logP Q.*/

#include "attacksECC.h"

using namespace Givaro;
using namespace std;

#define NumberOfPartitions_Pollard 16
#define pollardThreshold 10 //Pollard rho would try 10 times
#define noMatchThreshold 50
    
Integer& H(Integer& result,ecPoint& input, int L, ellipticCurveFq& E_Fq)
{
	Integer mod=(Integer)L;
	//Integer::random_lessthan (result,mod);
	ExtensionField::Element P,R=input.x,Q=input.y;
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
    if (Q.size()) 
	{
        E_Fq.field->Fp_X.assign(P, Q);
        E_Fq.field->Fp_X.setdegree(P);
        for(unsigned long l=0;l<P.size();++l) 
        {
    		result+=P[l]*(l+1);//if P[l] loses precision, that's not a problem
            result=result%mod;
        }                    
    }

    if(result<0)
    	result+=mod;
    return result;
}
//compute result such that Q=result*P
void pollardRho(Integer& result, ecPoint& P, ecPoint& Q,Integer n,ellipticCurveFq& E_Fq,int L=0)
{
    int count=0, noMatchCount=0;
	if(Q.identity)
	{
		result=0;
		return;
	}
	if(P==Q||n==2)
	{
		result=1;
		return;	
	}
	
	//std::cout<<"order of P0 is"<<n<<std::endl;
	if(L==0)
	{
		std::cout<<"No. of partitions L?";
		std::cin>>L;
	}
	Integer a[L],b[L];
	Integer c1,c2,d1,d2;
	//a=(Integer*)malloc(sizeof(Integer)*L);
	//b=(Integer*)malloc(sizeof(Integer)*L);
	
	ecPoint R[L],P1,P2,tP;
	Integer::seeding();
		
	for(int i=0; i<L; i++)
	{
		Integer::random_lessthan (a[i],n);
		Integer::random_lessthan (b[i],n);
		E_Fq.compute(R[i],a[i],P,b[i],Q,n);//R[i]=a[i]*P+b[i]*Q
	}

	Random:
    count++;
    if(count>pollardThreshold)
    {
        cout<<"pollard rho method failed to compute, try giving your own initial number c1 and d1"<<endl;
        exit(0);
    }
	Integer::random_lessthan (c1,n);
	Integer::random_lessthan (d1,n);
	
	E_Fq.compute(P1,c1,P,d1,Q,n);
	P2=P1;
	c2=c1;
	d2=d1;
	Integer j;
	do
	{
        noMatchCount++;
        H(j,P1, L, E_Fq);
		//std::cout<<"j: "<<j;
		E_Fq.add(tP,P1,R[j]);
		P1=tP;
		c1=(c1+a[j])%n;
		d1=(d1+b[j])%n;
		for(int i=0;i<2;i++)
		{	
			H(j,P2, L, E_Fq);
			//std::cout<<"j: "<<j;
			E_Fq.add(tP,P2,R[j]);
			P2=tP;
			c2=(c2+a[j])%n;
			d2=(d2+b[j])%n;	
		}
        if (noMatchCount>noMatchThreshold)
        {
            cout<<"pollard rho method failed to compute, try changing partition function"<<endl;
            exit(0);
        }
	}while(!(P2==P1));
	
	if(d1==d2)
	{
		//cout<<"trying pollard rho again";
		/*result=-1;*/
		
		goto Random;
	}
	else
	{
		Integer numr=c1-c2, denr=d2-d1;
		Integer d=gcd(denr,numr);
		denr=denr/d; numr=numr/d;
		d=gcd(denr, n);
		if(d==1)
		{
			invin (denr,n);
			result= ((numr*denr)%n);
			if(result<0)
				result+=n;
			//std::cout<<d<<","<<c1<<","<<c2<<","<<d1<<","<<d2<<":"<<result;
		}
		else
		{
		//cout<<"trying pollard rho again";
		goto Random;
		}
		
	}/**/
	
}
/*assumption: elements of p are pairwise relatively prime*/
void CRT(Integer& result,std::vector<Integer>& a,std::vector<Integer>& p)
{
    Integer length=p.size();
    Integer n=1;
    for(Integer i=0;i<length;i++)
        n*=p[i];
    std::vector<Integer> m(length),c(length);
    for(Integer i=0;i<length;i++)
    {
        m[i]=n/p[i];
        inv(c[i],m[i],p[i]);//mi^-1 mod p[i]
        c[i]*=m[i];
        result=(result+a[i]*c[i])%n;
    }
}

void pohligHellman(Integer& result, ecPoint& P, ecPoint& Q,Integer n,ellipticCurveFq& E_Fq)
{
    IntFactorDom<> IP;
    
    std::vector<Integer> Lf;
    std::vector<unsigned long> Le;
    IP.set(Lf,Le,n,0);
    Integer length=Lf.size();
    cout<<length;
    std::vector<Integer> res(length);
	ecPoint P0,_Q,temp1,temp2,S;
	Integer z, power,f;
    for(Integer i=0;i<length;i++)
    {
    	S=Q;
    	power=1;
        //cout<<Lf[i]<<"^"<<Le[i]<<"*";
        f=n/Lf[i];

        E_Fq.scalarMultiply(P0,P,f,n);
        E_Fq.scalarMultiply(_Q,S,f,n);//S is multiple of P, ord(P)=n
        //std::cout<<f<<"*"; E_Fq.show(S); std::cout<<"=";
        //E_Fq.show(P0); E_Fq.show(_Q);
        pollardRho(z,P0,_Q,Lf[i],E_Fq,NumberOfPartitions_Pollard);//ord(P0)=Lf[i]
        res[i]=z;
        //cout<<"res"<<i<<":"<<z<<std::endl;
        for(unsigned long j=1; j<Le[i];j++)
        {
            f=f/Lf[i];
            E_Fq.scalarMultiply(temp1,P,z*power,n);//ord(P)=n
            E_Fq.inv(temp2,temp1);
            E_Fq.add(temp1,S,temp2);
            S=temp1;
            E_Fq.scalarMultiply(_Q,S,f,n);//S is multiple of P, ord(P)=n
            pollardRho(z,P0,_Q,Lf[i],E_Fq,NumberOfPartitions_Pollard);//ord(P0)=Lf[i]
            //E_Fq.show(_Q);
            //std::cout<<"for i and j, z is"<<i<<","<<j<<","<<z<<std::endl;
            if(z==-1)
            {
            	return;
            }
            power=power*Lf[i];
            res[i]+=z*power;
        }
        Lf[i]=Lf[i]*power;  
    }
    /*cout<<"a%p"<<endl;
    for(Integer i=0;i<length;i++)
    {
    	cout<<res[i]<<":"<<Lf[i]<<endl;
    }*/
    CRT(result,res,Lf);//should send Lf[i]^Le[i]
    cout<<"Pohlig computes discrete logarithm LogpQ as "<<result<<endl;
}
