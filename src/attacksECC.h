#ifndef EC_H
#define EC_H
	#include "ellipticCurve.h"
#endif

#ifndef EC_GIVAROINTEGER_H
#define EC_GIVAROINTEGER_H
	#include "gmp++/gmp++.h"
	#include <givaro/givinit.h>
	#include <givaro/givintfactor.h>
	#include <givaro/givtimer.h>
#endif

using namespace Givaro;

/*computes result such that Q=result*P. 
Assumption n is prime, n=order(P) in E(Fq) represented by E_Fq, L is number of partition functions*/
void pollardRho(Integer& result, ecPoint& P, ecPoint& Q,Integer n,ellipticCurveFq& E_Fq,int L);
/*computes result such that Q=result*P. 
n=order(P) in E(Fq) represented by E_Fq. 
It uses function pollardRho to compute each smaller instance of ECDLP*/
void pohligHellman(Integer& result, ecPoint& P, ecPoint& Q,Integer n,ellipticCurveFq& E_Fq);

/*Auxillary Functions*/

/*result=H(input) where 0=<result<L. used by function pollardRho to compute the partition*/
Integer& H(Integer& result,ecPoint& input, int L, ellipticCurveFq& E_Fq);
/*Chinese Remainder theorem, used by function pohligHellman to compute result, where result=a[i]%p[i]*/
void CRT(Integer& result,std::vector<Integer>& a,std::vector<Integer>& p);



