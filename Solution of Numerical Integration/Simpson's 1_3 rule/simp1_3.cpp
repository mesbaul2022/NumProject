#include<bits/stdc++.h>

using namespace std;
double func(double x)
{
    return 2/(1+pow(x,3));
}
int main()
{
    cout<<"Enter the value of a b n"<<endl;
    double a,b,n;
    cin>>a>>b>>n;
    double h=(b-a)/n;
    double sum=func(a)+func(b);
    for(int i=1;i<n;i++)
    {
        double x=a+(i*h);
        if(i%2==0)
        {
            sum+=2*func(x);
        }
        else
        {
            sum+=4*func(x);
        }
    }
    cout<<"result"<<endl;
    cout<<(sum*h)/3<<endl;

}