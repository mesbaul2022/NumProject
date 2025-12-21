#include<bits/stdc++.h>
using namespace std;
int fact(int n)
{
    if(n==1 || n==0)
    {
        return 1;
    }
    return n*fact(n-1);
}

int main()
{
    int n;
    cout<<"Enter the n :";
    cin>>n;
    vector<double>X(n);
    vector<vector<double>>Y(n,vector<double>(n,0.0));
    for(int i=0;i<n;i++)
    {
        cin>>X[i];
    }
    cout<<"Enter the Y "<<endl;
     for(int i=0;i<n;i++)
    {
        cin>>Y[i][0];
    }

    cout<<"TABLE"<<endl;

    for(int j=1;j<n;j++)
    {
        for(int i=n-1;i>=j;i--)
        {
           Y[i][j]=Y[i][j-1]-Y[i-1][j-1];
        }
    }

    for(int i=0;i<n;i++)
    {
        cout<<X[i]<<"\t";
        for(int k=0;k<=i;k++)
        {
           cout<<Y[i][k]<<"\t";
        }
        cout<<endl;
    }

    double h=X[1]-X[0];
    double value;
    cout<<"Enter the value : ";
    cin>>value;
    double v=(value - X[n-1])/h;
    double result=Y[n-1][0];
    double t=1.0;
    for(int i=1;i<n;i++)
    {
           t*=(v+(i-1));
           result+=(t*Y[n-1][i])/fact(i);
    }
    cout<<"The final result "<<result<<endl;
    return 0;
}