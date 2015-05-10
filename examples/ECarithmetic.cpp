#include "../src/ellipticCurve.h"

using namespace Givaro;
using namespace std;

int main()
{
    ellipticCurveFq E_Fq;
    E_Fq.show();
    ExtensionField::Element x1,y1,x2,y2,x3,y3;
    cout<<"Write x1: ";
    E_Fq.field->readElement(x1,true);
    cout<<"Write y1: ";
    E_Fq.field->readElement(y1,true);
    cout<<"Write x2: ";
    E_Fq.field->readElement(x2,true);
    cout<<"Write y2: ";
    E_Fq.field->readElement(y2,true);
    ecPoint P(x1,y1);
    ecPoint Q(x2,y2);
    ecPoint R;
    E_Fq.Double(R,P);//R=2*P
    cout<<endl;
    cout<<"2*"; E_Fq.show(P); cout<<" = "; E_Fq.show(R); cout<<endl;
    E_Fq.add(R,P,Q);//R=P+Q
    E_Fq.show(P); cout<<" + "; E_Fq.show(Q); cout<<" = ";  E_Fq.show(R); cout<<endl;
    E_Fq.scalarMultiply(R,P,(Integer)6,-1);//R=6*P, order of P is not required
    cout<<"6*"; E_Fq.show(P); cout<<" = "; E_Fq.show(R); cout<<endl;
    cout<<"Input (x,y) to verify if it is a Point on E(Fq): ";
    E_Fq.field->readElement(x3, true);
    E_Fq.field->readElement(y3, true);
    ecPoint T(x3,y3); 
    if(E_Fq.verifyPoint(T))
    {
        cout<<"T is a point.";
    } 
    else
    {
        cout<<"T is not a point.";
    }  
    cout<<endl;
    E_Fq.inv(R,P);
    cout<<"Inverse of "; E_Fq.show(P); cout<<" = "; E_Fq.show(R);
    cout<<endl;   
    return 0;
}
