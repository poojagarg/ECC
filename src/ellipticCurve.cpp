#ifndef EC_H
#define EC_H
    #include "ellipticCurve.h"
#endif

using namespace Givaro;
using namespace std;

ExtensionField::ExtensionField(Integer p, Integer m)
{
    this->p=p;
    this->m=m;
    Fp=primeField((primeField::Residu_t)p);
    Fp_X=Poly1Dom< primeField, Dense >( Fp, Indeter("X") );
    if(m==1)
        Fp_X.assign(irred,(Element)p);
    else
    {
        
        cout<<"Give an irreducible polynomial of degree "<<m<<": ";
        Fp_X.read(cin,irred);
        cout<<endl;
    }
    primeField::Element tmp;
    Fp_X.assign(zero,Fp_X.zero);
    Fp_X.assign(one,Fp_X.one);
    Fp_X.assign(mOne,Fp_X.mOne);
    
}
    
ExtensionField ExtensionField::operator=(const ExtensionField& E)
{
    p=E.p;
    m=E.m;
    Fp=E.Fp;
    Fp_X=E.Fp_X;
    irred=E.irred;
    zero=E.zero;
    one=E.one;
    mOne=E.mOne;
    return *this;
}
//Return R=(A+B)%p
ExtensionField::Element& ExtensionField::add(Element& R, const Element& A, const Element& B) const
{
    Fp_X.add(R,A,B);
    return R;
}
ExtensionField::Element& ExtensionField::addin(Element& R, const Element& A) const
{
    Fp_X.addin(R,A);
    return R;
}
//Return R=(A-B)%p
ExtensionField::Element& ExtensionField::sub(Element& R, const Element& A, const Element& B) const
{
    Fp_X.sub(R,A,B);
    return R;
}
//Return R=-B
ExtensionField::Element& ExtensionField::neg(Element& R, const Element& A) const
{
    Fp_X.neg(R,A);
    return R;
}
//Return R=A*B%irred
ExtensionField::Element& ExtensionField::mul(Element& R, const Element& A, const Element& B) const
{
    if(m>1)
    {
        Element tmp,Q;
        Fp_X.mul(tmp,A,B);
        Fp_X.divmod(Q,R,tmp,irred);
    }
    else
        Fp_X.mul(R,A,B);
    return R;
}
//Return R=A*A%irred
ExtensionField::Element& ExtensionField::sqr(Element& R, const Element& A) const
{
    if(m>1)
    {
        Element tmp,Q;
        Fp_X.sqr(tmp,A);
        Fp_X.divmod(Q,R,tmp,irred);
    }
    else
        Fp_X.sqr(R,A);
    return R;
}
//Return inverse of (A) with irred polynomial, if gcd(A,irred)!=1, then it return 0
ExtensionField::Element& ExtensionField::inv(Element& I, const Element& A) const
{
    Element tmp,D;
    D=Fp_X.gcd(D,I,tmp,A,irred);
    if(Fp_X.isOne(D))
        return I;
    Fp_X.assign(I,Fp_X.zero);
    return I;

}
ExtensionField::Element& ExtensionField::additiveInv(Element& I, const Element& A) const
{
    Fp_X.sub(I,Fp_X.zero,A);
    return I;

}

//Return Q, such that Q=A/B=A*(inverse of B),
ExtensionField::Element& ExtensionField::div(Element& Q, const Element& A, const Element& B) const
{
    Element I;
    I=inv(I,B);
    mul(Q,A,I);
    return Q;
}
ExtensionField::Element& ExtensionField::scalarMultiply(Element& R,const Element& A,Integer k) const
{
    k=k%p;
    if(m==1)
    {
        Element eK;
        primeField::Element tmp;
        Fp_X.assign(eK,Fp.init(tmp,k));
        Fp_X.mul(R,A,eK);
        return R;
    }
    Element tmp,acc=A;
    Fp_X.assign(R,zero);
    while(k>0)
    {
        if(k%2)
            Fp_X.addin(R,acc);
        acc=Fp_X.add(tmp,acc,acc); //acc=tmp=2*acc
        k=k/2;
    }
    return R;
}
void ExtensionField::readElement(Element& P, bool flag)
{
    Fp_X.read(cin,P,flag);
}
void ExtensionField::writeElement(Element& A)
{
    cout<<"( ";
    Fp_X.write(cout<<"",A); 
    cout<<" )";  
}
bool ExtensionField::isElement(const Element& A)
{
	Degree d=Fp_X.degree(A);
	if(d._deg<m)
	    return true;
	return false;
}

ecPoint::ecPoint(bool b)
{
    identity=b;
}
ecPoint::ecPoint(Element x1, Element y1)
: identity(false) {
    x=x1;
    y=y1;
}
ecPoint ecPoint::operator=( const ecPoint& P)
{
    if(P.identity)
    {
        this->identity=true;
        return *this;
    }
    
    this->identity=P.identity;
    this->x=P.x;
    this->y=P.y;
    return *this;
}
bool ecPoint::operator==(const ecPoint& t) const
{return (identity && t.identity) || (!identity && !t.identity && x==t.x && y==t.y);}
bool ecPoint::operator< (const ecPoint& t) const
{return identity ? !t.identity : (!t.identity && (x<t.x || (x==t.x && y<t.y)));}



ellipticCurve::ellipticCurve()
{
    Integer p, m;
    cin>>p>>m;
    Kptr=new ExtensionField(p,m);
    cout<<"Write A: ";
    Kptr->readElement( A);
    cout<<"Write B: ";
    Kptr->readElement( B);
    cout<<"type? ";
    cin>>type;
    if(type==2)
    {
        cout<<"Write C: ";
        Kptr->readElement( C);
    }
    cout<<endl;
}

ellipticCurve::~ellipticCurve()
{
    free(Kptr);
}

ellipticCurve ellipticCurve::operator=( const ellipticCurve& F)
{
    Kptr=F.Kptr;
    A=F.A;
    B=F.B;
    C=F.C;
    type=F.type;
    return *this;
}
void ellipticCurve::print()
{
    cout<<"A:";
    Kptr->writeElement(A);
    cout<<endl;
    cout<<"B:";
    Kptr->writeElement(B);
    cout<<endl;
    cout<<"C:";
    Kptr->writeElement(C);
    cout<<endl;
    cout<<"type"<<type<<endl;
}

ellipticCurveFq::ellipticCurveFq()
{
    ec=new ellipticCurve();
    cout<<"degree d? ";
    cin>>d;
    if(d!=1)
        field=new ExtensionField(ec->Kptr->p,(ec->Kptr->m)*d);
    else
        field=ec->Kptr;
    identity.identity=true;
    this->d=d;
}
ellipticCurveFq::ellipticCurveFq(ellipticCurve* e)
{
    ec=e;
    cout<<"degree d? ";
    cin>>d;
    if(d!=1)
        field=new ExtensionField(ec->Kptr->p,(ec->Kptr->m)*d);
    else
        field=ec->Kptr;
    identity.identity=true;
    this->d=d;
}
ellipticCurveFq::~ellipticCurveFq()
{  
    free(ec);
    free(field);
}
const ellipticCurveFq::Point& ellipticCurveFq::inv(Point& Q, const Point &P) 
{
    if(P.identity)
    {
        Q.identity=true;
        return Q;
    }
    Q.identity=false;
    Point T;
    switch(ec->type)
    {
        case 0:
            Q.x=P.x;
            field->additiveInv(Q.y,P.y); //Q.y+P.y=0
            return Q;
            
        case 1:
            Q.x=P.x;
            field->add(Q.y,P.y,P.x);
            //Q.y+=Q.x;
            return Q;
            
        case 2:
            Q.x=P.x;
            field->add(Q.y,P.y,ec->C);

            //Q.y+=ec->C;
            return Q;
        default:
            return Q;
            
    }
}
bool ellipticCurveFq::isInv(const Point& Q, const Point &P) //is Q+P=point at inifinity?
{
    Point R;
    inv(R,P);
    //cout<<"At inverse check"; show(R);
    if(R==Q)
        return true;
    return false;
}

ellipticCurveFq::Point& ellipticCurveFq::Double(Point &R,Point &P) 
{
    if (P.identity||isInv(P,P)) 
    {
        R.identity=true;
        return R;
    }

    R.identity=false;
    fieldElement x,y,x12;
    fieldElement slopeSquare, slope;
    field->sqr(x12,P.x);
    fieldElement _3x12,_2y,_3x12pA,_2x,_x,slope_x;
    fieldElement xpy,Bdx12,slopex,slopex_x;
    fieldElement x12pA,ypC,xpx,slopexpx;
    switch(ec->type)
    {
        case 0:
            
            field->scalarMultiply(_3x12,x12,3);
            field->scalarMultiply(_2y,P.y,2);
            field->add(_3x12pA,_3x12,ec->A);
            field->div(slope,_3x12pA,_2y);
            field->sqr(slopeSquare,slope);
            field->scalarMultiply(_2x,P.x,2);
            field->sub(x,slopeSquare,_2x);
            field->sub(_x,P.x,x);
            field->sub(y,field->mul(slope_x,slope,_x),P.y);
        
            R.x=x;
            R.y=y;
            return R;

        case 1:

            field->add(xpy,x12,P.y);
            field->div(slope,xpy,P.x);
            field->div(Bdx12,ec->B,x12);
            field->add(x,x12,Bdx12);
            field->mul(slopex,slope,x);
            field->add(slopex_x,slopex,x);
            field->add(y,x12,slopex_x);
            R.x=x;
            R.y=y;
            return R;
        case 2:
            
            field->add(x12pA,x12,ec->A);
            field->div(slope,x12pA,ec->C);
            field->sqr(slopeSquare,slope);
            x=slopeSquare;
            field->add(ypC,P.y,ec->C);
            field->add(xpx,P.x,x);
            field->mul(slopexpx,xpx,slope);
            field->add(y,ypC,slopexpx);
            R.x=x;
            R.y=y;
            return R;
    } 
    return R;
}
ellipticCurveFq::Point& ellipticCurveFq::add(Point &R,Point &P, Point &Q) 
{
    if(P.identity&&Q.identity||isInv(P,Q))
    {
        R.identity=true;
        return R;
    }
    R.identity=false;
    if (P.identity) 
    {
        R.x=Q.x;
        R.y=Q.y;
        return R;
    }
    if (Q.identity) 
    {
        R.x=P.x;
        R.y=P.y;
        return R;
    }
    if(P==Q)
        return Double(R,P);
    
    fieldElement x,y,x12;
    fieldElement slopeSquare, slope;
    fieldElement y2m1,x2m1,x1p2,x1m3,slopex1m3,y1pC,x1p2pA,x1p2pApS,slopex1m3px3;
    field->sub(y2m1,Q.y,P.y);
    field->sub(x2m1,Q.x,P.x);
    field->div(slope,y2m1,x2m1);
    field->sqr(slopeSquare,slope);
    field->add(x1p2,P.x,Q.x);
    switch(ec->type)
    {
        case 0:
            field->sub(x,slopeSquare,x1p2);
            field->sub(x1m3,P.x,x);
            field->mul(slopex1m3,slope,x1m3);
            field->sub(y,slopex1m3,P.y);
            R.x=x;
            R.y=y;
            return R;
        case 1:
            field->add(x1p2pA,x1p2,ec->A);
            field->add(x1p2pApS,x1p2pA,slope);
            field->add(x,x1p2pApS,slopeSquare);
            field->add(x1m3,P.x,x);
            field->mul(slopex1m3,slope,x1m3);
            field->add(slopex1m3px3,slopex1m3,x);
            field->add(y,slopex1m3px3,P.y);
            R.x=x;
            R.y=y;
            return R;
        case 2:
            field->add(x,slopeSquare,x1p2);
            field->add(x1m3,P.x,x);
            field->mul(slopex1m3,slope,x1m3);
            field->add(y1pC,P.y,ec->C);
            field->add(y,y1pC,slopex1m3);
            R.x=x;
            R.y=y;
            return R;
    } 
    return R;
}
void ellipticCurveFq::show(Point& P)
{
    if(P.identity)
    {
        cout<<"o";
        return;
    }

    cout<<"(";
    field->writeElement(P.x);
    cout<<", ";  
    field->writeElement(P.y);  
    cout<<")";
}
//R=kP
ellipticCurveFq::Point& ellipticCurveFq::scalarMultiply(Point&R, Point& P, Integer k, Integer order)//order of P 
{
    if(P.identity)
    {
        R.identity=true;
        return R;
    }
    Point temp(false);
    R.identity=true;
    Point acc(P);
    if(order!=-1)//if order is known
        k=k%order;//better: k=k%order(P)
    while(k>0)
    {
        if(k%2)
        {
            /*if(k==1)
            {
                show(R);
                show(acc);
            }*/
            add(temp,R,acc);
            R=temp;
        }
        Double(temp,acc);    
        acc=temp;
        k=k/2;
    }
    return R;
}
bool ellipticCurveFq::verifyPoint(const Point &P) const
{
    if(P.identity)
        return true;
    fieldElement rhs,x2,y2,x2pAx,x2pA,xpA,xpAx2,lhs,xy,cy;
    field->sqr(x2,P.x);
    field->sqr(y2,P.y);
    
    switch (ec->type) 
    {
        case 0:
            field->add(x2pA,x2,ec->A);
            field->mul(x2pAx,x2pA,P.x);
            field->add(rhs,x2pAx,ec->B);
            /*field->writeElement(rhs)<<endl;  
            field->writeElement(y2)<<endl;  
            */
            if(y2==rhs)
                return true;
            else 
                return false;
        case 1:
            field->add(xpA,P.x,ec->A);
            field->mul(xpAx2,xpA,x2);
            field->add(rhs,xpAx2,ec->B);
            field->mul(xy,P.x,P.y);
            field->add(lhs,y2,xy);
            if(lhs==rhs)
                return true;
            else 
                return false;
        case 2:
            field->add(x2pA,x2,ec->A);
            field->mul(x2pAx,x2pA,P.x);
            field->add(rhs,x2pAx,ec->B);
            field->mul(cy,ec->C,P.y);
            field->add(lhs,y2,cy);
            if(lhs==rhs)
                return true;
            else 
                return false;
        default:
            return false;
    }
}
/*type 0: E/K, char(K)!=2: y2 = x3+ax+b,
 type 1: non-supersingular E/F2m: y2 + xy = x3+ax2+b,
 tyep 2: supersingular E/F2m: y2 + cy = x3 + ax + b*/
void ellipticCurveFq::show()
{
    cout<<"Elliptic Curve Defined by ";
    switch(ec->type)
    {
        case 0: 
            cout<<"y^2 = x^3 + "; 
            field->writeElement(ec->A);
            cout<<"x + ";
            field->writeElement(ec->B);
            cout<<endl;
            break;
        case 1:
            cout<<"y^2 + xy = x^3 + "; 
            field->writeElement(ec->A);
            cout<<"x^2 + ";
            field->writeElement(ec->B);
            cout<<endl;
            break;
        case 2:
            cout<<"y^2 + "; 
            field->writeElement(ec->C);
            cout<<"y = x^3 + ";
            field->writeElement(ec->A);
            cout<<"x + ";
            field->writeElement(ec->B);
            cout<<endl;
            break;
    }
    cout<<" over finite field in X of size "<<field->p<<"^"<<field->m;
    cout<<" with irreducible polynomial ";
    field->writeElement(field->irred);
    cout<<endl;
}
void ellipticCurveFq::compute(ecPoint& R,Integer a, ecPoint& P,Integer b,ecPoint& Q,Integer n)
{
    ecPoint t1,t2;
    scalarMultiply(t1,P,a,n);
    scalarMultiply(t2,Q,b,n);
    add(R,t1,t2);
}
