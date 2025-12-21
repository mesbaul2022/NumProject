#include<bits/stdc++.h>
using namespace std;
double f(double x,double y){
  return (x-y)/2.0;
}
int main(){
cout<<" x0 and y0 :";
  double x0,y0;
  cin>>x0>>y0;
  double h;
  cout<<" h:";
   cin>>h;
   double xn;
   cout<<"xn :";
   cin>>xn;
   double interval=(xn-x0)/h;
   double yn;
   for(double i=1;i<=interval;i++){
       double k1=h*f(x0,y0);
       double k2=h*f(x0+(h/2.0),y0+(k1/2.0));
       double k3=h*f(x0+(h/2.0),y0+(k2/2.0));
       double k4=h*f(x0+h,y0+k3);
       double k=(k1+2*k2+2*k3+k4)/6;
       yn=y0+k;
      cout  <<"x0= "<< x0<<"  y0= "<< y0<<"  yn="<<yn<<endl;
       y0=yn;
       x0+=h;

   }
   cout<<"y("<<xn<<")= "<<yn<<endl;
return 0;
}
