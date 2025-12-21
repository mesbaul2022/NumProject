# Numerical Lab Project (C++)

Numerical Methods Laboratory implementations in **C++**: solution of linear equations, non-linear root finding, interpolation, numerical integration & differentiation, least squares regression, and Runge–Kutta (ODE).

---

##Project Overview:
- [Overview](overview)
---

## Table of Contents

- [Solution of Linear Equations](#solution-of-linear-equations)
  - [Gauss Elimination Method](#gauss-elimination)
    - [Theory](#gauss-elimination-theory)
    - [Code](#gauss-elimination-code)
    - [Input](#gauss-elimination-input)
    - [Output](#gauss-elimination-output)
  - [Gauss Jordan Elimination Method](#gauss-jordan)
    - [Theory](#gauss-jordan-theory)
    - [Code](#gauss-jordan-code)
    - [Input](#gauss-jordan-input)
    - [Output](#gauss-jordan-output)
  - [LU Decomposition Method](#lu-decomposition)
    - [Theory](#lu-decomposition-theory)
    - [Code](#lu-decomposition-code)
    - [Input](#lu-decomposition-input)
    - [Output](#lu-decomposition-output)
  - [Matrix Inversion](#matrix-inversion)
    - [Theory](#matrix-inversion-theory)
    - [Code](#matrix-inversion-code)
    - [Input](#matrix-inversion-input)
    - [Output](#matrix-inversion-output)

- [Solution of Non-Linear Equations](#solution-of-non-linear-equations)
  - [Bisection Method](#bisection)
    - [Theory](#bisection-theory)
    - [Code](#bisection-code)
    - [Input](#bisection-input)
    - [Output](#bisection-output)
  - [False Position Method](#false-position)
    - [Theory](#false-position-theory)
    - [Code](#false-position-code)
    - [Input](#false-position-input)
    - [Output](#false-position-output)
  - [Secant Method](#secant)
    - [Theory](#secant-theory)
    - [Code](#secant-code)
    - [Input](#secant-input)
    - [Output](#secant-output)
  - [Newton Raphson Method](#newton-raphson)
    - [Theory](#newton-raphson-theory)
    - [Code](#newton-raphson-code)
    - [Input](#newton-raphson-input)
    - [Output](#newton-raphson-output)

- [Solution of Interpolation](#interpolation)
  - [Newton's Forward Interpolation Method](#newtons-forward-interpolation-method)
    - [Theory](#newtons-forward-interpolation-theory)
    - [Code](#newtons-forward-interpolation-code)
    - [Input](#newtons-forward-interpolation-input)
    - [Output](#newtons-forward-interpolation-output)
  - [Newton's Backward Interpolation Method](#newtons-backward-interpolation-method)
    - [Theory](#newtons-backward-interpolation-theory)
    - [Code](#newtons-backward-interpolation-code)
    - [Input](#newtons-backward-interpolation-input)
    - [Output](#newtons-backward-interpolation-output)
  - [Divided Difference Method](#divided-difference-method)
    - [Theory](#divided-difference-theory)
    - [Code](#divided-difference-code)
    - [Input](#divided-difference-input)
    - [Output](#divided-difference-output)

- [Solution of Curve Fitting Model](#solution-of-curve-fitting-model)
  - [Least Square Regression Method For Linear Equations](#least-square-regression-method-for-linear-equations-method)
    - [Theory](#least-square-regression-method-for-linear-equations-theory)
    - [Code](#least-square-regression-method-for-linear-equations-code)
    - [Input](#least-square-regression-method-for-linear-equations-input)
    - [Output](#least-square-regression-method-for-linear-equations-output)
  - [Least Square Regression Method For Transcendental Equations](#least-square-regression-method-for-transcendental-equations)
    - [Theory](#least-square-regression-method-for-transcendental-equations-theory)
    - [Code](#least-square-regression-method-for-transcendental-equations-code)
    - [Input](#least-square-regression-method-for-transcendental-equations-input)
    - [Output](#least-square-regression-method-for-transcendental-equations-output)
  - [Least Square Regression Method For Polynomial Equations](#least-square-regression-method-for-polynomial-equations)
    - [Theory](#least-square-regression-method-for-polynomial-equations-theory)
    - [Code](#least-square-regression-method-for-polynomial-equations-code)
    - [Input](#least-square-regression-method-for-polynomial-equations-input)
    - [Output](#least-square-regression-method-for-polynomial-equations-output)

- [Solution of Differential Equations](#solution-of-differential-equations)
  - [Runge Kutta Method](#runge-kutta-method)
    - [Theory](#runge-kutta-theory)
    - [Code](#runge-kutta-code)
    - [Input](#runge-kutta-input)
    - [Output](#runge-kutta-output)

- [Numerical Differentiation](#numerical-differentiation)
  - [Numerical Differentiation by Forward Interpolation Method](#numerical-differentiation-by-forward-interpolation-method)
    - [Theory](#numerical-differentiation-by-forward-interpolation-theory)
    - [Code](#numerical-differentiation-by-forward-interpolation-code)
    - [Input](#numerical-differentiation-by-forward-interpolation-input)
    - [Output](#numerical-differentiation-by-forward-interpolation-output)
  - [Numerical Differentiation by Backward Interpolation Method](#numerical-differentiation-by-backward-interpolation-method)
    - [Theory](#numerical-differentiation-by-backward-interpolation-theory)
    - [Code](#numerical-differentiation-by-backward-interpolation-code)
    - [Input](#numerical-differentiation-by-backward-interpolation-input)
    - [Output](#numerical-differentiation-by-backward-interpolation-output)

- [Solution of Numerical Integrations](#solution-of-numerical-integrations)
  - [Simpson's One-Third Rule](#simpson-13)
    - [Theory](#simpsons-one-third-rule-theory)
    - [Code](#simpsons-one-third-rule-code)
    - [Input](#simpsons-one-third-rule-input)
    - [Output](#simpsons-one-third-rule-output)
  - [Simpson's Three-Eighths Rule](#simpson-38)
    - [Theory](#simpsons-three-eighths-rule-theory)
    - [Code](#simpsons-three-eighths-rule-code)
    - [Input](#simpsons-three-eighths-rule-input)
    - [Output](#simpsons-three-eighths-rule-output)
  


---

## Solution of Linear Equations

### Gauss Elimination

**Goal:** Solve \(Ax=b\) using forward elimination + back substitution.  
**File:** `Linear/gauss_elimination.cpp`

**Input format**
- `n`
- `A` matrix (`n x n`)
- `b` vector (`n`)

```
#include <bits/stdc++.h>
using namespace std;

static const double EPS = 1e-12;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;
    vector<vector<double>> a(n, vector<double>(n+1));
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++) cin >> a[i][j];
    }
    for(int i=0;i<n;i++) cin >> a[i][n];

    // Partial pivoting + forward elimination
    for(int col=0; col<n; col++){
        int pivot = col;
        for(int row=col+1; row<n; row++)
            if (fabs(a[row][col]) > fabs(a[pivot][col])) pivot = row;

        if (fabs(a[pivot][col]) < EPS) {
            cout << "No unique solution (singular matrix).\n";
            return 0;
        }
        swap(a[pivot], a[col]);

        for(int row=col+1; row<n; row++){
            double factor = a[row][col] / a[col][col];
            for(int k=col; k<=n; k++){
                a[row][k] -= factor * a[col][k];
            }
        }
    }

    // Back substitution
    vector<double> x(n);
    for(int i=n-1; i>=0; i--){
        double sum = a[i][n];
        for(int j=i+1; j<n; j++) sum -= a[i][j]*x[j];
        x[i] = sum / a[i][i];
    }

    cout.setf(std::ios::fixed); cout<<setprecision(6);
    for(int i=0;i<n;i++) cout << "x["<<i<<"] = " << x[i] << "\n";
    return 0;
}
```

---

### Gauss Jordan

**Goal:** Solve \(Ax=b\) by reducing \([A|b]\) to RREF.  
**File:** `Linear/gauss_jordan.cpp`

```
#include <bits/stdc++.h>
using namespace std;

static const double EPS = 1e-12;

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;
    vector<vector<double>> a(n, vector<double>(n+1));
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++) cin >> a[i][j];
    }
    for(int i=0;i<n;i++) cin >> a[i][n];

    for(int col=0; col<n; col++){
        int pivot = col;
        for(int row=col; row<n; row++)
            if (fabs(a[row][col]) > fabs(a[pivot][col])) pivot = row;

        if (fabs(a[pivot][col]) < EPS){
            cout << "No unique solution (singular matrix).\n";
            return 0;
        }
        swap(a[pivot], a[col]);

        double div = a[col][col];
        for(int k=col; k<=n; k++) a[col][k] /= div;

        for(int row=0; row<n; row++){
            if(row==col) continue;
            double factor = a[row][col];
            for(int k=col; k<=n; k++){
                a[row][k] -= factor * a[col][k];
            }
        }
    }

    cout.setf(std::ios::fixed); cout<<setprecision(6);
    for(int i=0;i<n;i++) cout << "x["<<i<<"] = " << a[i][n] << "\n";
    return 0;
}
```

---

### LU Decomposition

**Goal:** Factorize \(A = LU\) (Doolittle) and solve \(Ax=b\).  
**File:** `Linear/lu_decomposition.cpp`

```
#include <bits/stdc++.h>
using namespace std;

static const double EPS = 1e-12;

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;
    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) cin >> A[i][j];
    for(int i=0;i<n;i++) cin >> b[i];

    vector<vector<double>> L(n, vector<double>(n,0.0)), U(n, vector<double>(n,0.0));
    for(int i=0;i<n;i++) L[i][i] = 1.0;

    for(int i=0;i<n;i++){
        for(int j=i;j<n;j++){
            double sum = 0;
            for(int k=0;k<i;k++) sum += L[i][k]*U[k][j];
            U[i][j] = A[i][j] - sum;
        }
        if (fabs(U[i][i]) < EPS){
            cout << "LU failed (zero pivot). Try pivoting-based LU.\n";
            return 0;
        }
        for(int j=i+1;j<n;j++){
            double sum = 0;
            for(int k=0;k<i;k++) sum += L[j][k]*U[k][i];
            L[j][i] = (A[j][i] - sum) / U[i][i];
        }
    }

    // Forward substitution: Ly=b
    vector<double> y(n);
    for(int i=0;i<n;i++){
        double sum = b[i];
        for(int k=0;k<i;k++) sum -= L[i][k]*y[k];
        y[i] = sum; // since L[i][i]=1
    }

    // Back substitution: Ux=y
    vector<double> x(n);
    for(int i=n-1;i>=0;i--){
        double sum = y[i];
        for(int k=i+1;k<n;k++) sum -= U[i][k]*x[k];
        x[i] = sum / U[i][i];
    }

    cout.setf(std::ios::fixed); cout<<setprecision(6);
    for(int i=0;i<n;i++) cout << "x["<<i<<"] = " << x[i] << "\n";
    return 0;
}
```

---

### Matrix Inversion

**Goal:** Compute \(A^{-1}\) using Gauss-Jordan on \([A|I]\).  
**File:** `Linear/matrix_inversion.cpp`

```
#include <bits/stdc++.h>
using namespace std;

static const double EPS = 1e-12;

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n; 
    cin >> n;
    vector<vector<double>> aug(n, vector<double>(2*n, 0.0));

    // Read A
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++) cin >> aug[i][j];
    }
    // Append Identity
    for(int i=0;i<n;i++) aug[i][n+i] = 1.0;

    for(int col=0; col<n; col++){
        int pivot = col;
        for(int row=col; row<n; row++)
            if (fabs(aug[row][col]) > fabs(aug[pivot][col])) pivot = row;

        if (fabs(aug[pivot][col]) < EPS){
            cout << "Singular matrix. Inverse does not exist.\n";
            return 0;
        }
        swap(aug[pivot], aug[col]);

        double div = aug[col][col];
        for(int k=0;k<2*n;k++) aug[col][k] /= div;

        for(int row=0; row<n; row++){
            if(row==col) continue;
            double factor = aug[row][col];
            for(int k=0;k<2*n;k++){
                aug[row][k] -= factor * aug[col][k];
            }
        }
    }

    cout.setf(std::ios::fixed); cout<<setprecision(6);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cout << aug[i][n+j] << (j+1==n?'\n':' ');
        }
    }
    return 0;
}
```

---

## Solution of Non-linear Equations

> Edit the `f(x)` (and `df(x)` if needed) inside each file.

### Bisection

**File:** `NonLinear/bisection.cpp`

```
#include<bits/stdc++.h>
using namespace std;
int degree;
vector<double>coeff(degree+1);
 double f(double x){
        double val=0;
     double d=degree;
     int n=coeff.size();
     for(int i=0;i<n;i++){
           //  cout<<degree<<endl;
        double t=coeff[i]*(pow(x,d));
        val+=t;
        d--;
     }

   return val;
 }
int main(){

    cout<<"please1 enter the degree of polynomial equation :";
    cin>>degree;
 coeff.resize(degree + 1);
    cout<<"please enter the coefficient of polynomial :";
    for(int i=0;i<degree+1;i++)cin>>coeff[i];
    double xmax=sqrt((coeff[1]/coeff[0])*(coeff[1]/coeff[0])-2*(coeff[2]/coeff[0]));
  double a,b,c,root=0;
  for(double i=-xmax;i<=xmax;i+=0.5){
        a=i;
      b=i+0.5;
        double fa=f(a),fb=f(b);
   // cout<<a<<" "<<b<<" "<<fa*fb<<endl;
    if(fa*fb<0){
       double it=0;
       root++;
   double e=0.0001;
  do{
        it++;
    //cout<<"ok"<<endl;
    c= (a+b)/2.0;
   // cout<<"a="<<a<<" f(a)="<<f(a)<<" b="<<b<<" f(b)="<<f(b)<<" c="<<c<<" f(c)="<<f(c)<<endl;
    if(fabs(c)<=e) break;
    if(f(c)*f(a)<0)b=c;
    else a=c;
  }while(fabs(f(c))> e && fabs(b-a)>e);
   cout<<"the root "<< root<<" is: "<<c<<endl;
    cout<<"search interval ["<<i<<","<<i+0.5<<"]"<<endl;
    cout<<"iteration is:"<<it<<endl<<endl;
    }
  }




}

```

---

### False Position

**File:** `NonLinear/false_position.cpp`

```
#include<bits/stdc++.h>
using namespace std;
int degree;
vector<double>coeff(degree+1);
 double f(double x){
        double val=0;
     double d=degree;
     int n=coeff.size();
     for(int i=0;i<n;i++){
           //  cout<<degree<<endl;
        double t=coeff[i]*(pow(x,d));
        val+=t;
        d--;
     }

   return val;
 }
int main(){

    cout<<"please enter the degree of polynomial equation :";
    cin>>degree;
 coeff.resize(degree + 1);
    cout<<"please enter the coefficient of polynomial :";
    for(int i=0;i<degree+1;i++)cin>>coeff[i];
    double xmax=sqrt((coeff[1]/coeff[0])*(coeff[1]/coeff[0])-2*(coeff[2]/coeff[0]));
  double a,b,c,root=0;
  for(double i=-xmax;i<=xmax;i+=0.5){
        a=i;
      b=i+0.5;
        double fa=f(a),fb=f(b);
   // cout<<a<<" "<<b<<" "<<fa*fb<<endl;
    if(fa*fb<0){
       double it=0;
       root++;
   double e=0.0001;
  do{
        it++;

     c=a-f(a)*((b-a)/(f(b)-f(a)));
   // cout<<"a="<<a<<" f(a)="<<f(a)<<" b="<<b<<" f(b)="<<f(b)<<" c="<<c<<" f(c)="<<f(c)<<endl;
    if(fabs(c)<=e) break;
    if(f(c)*f(a)<0)b=c;
    else a=c;
  }while(fabs(f(c))> e && fabs(b-a)>e);
   cout<<"the root "<< root<<" is: "<<c<<endl;
    cout<<"search interval ["<<i<<","<<i+0.5<<"]"<<endl;
    cout<<"iteration is:"<<it<<endl<<endl;
    }
  }




}



```

---

### Secant

**File:** `NonLinear/secant.cpp`

```
#include<bits/stdc++.h>
using namespace std;
vector<double>coeff;
 void print(vector<double>coeff){
     int power=coeff.size()-1;
     bool m=false;
     for(int i=0;i<coeff.size();i++){
        if(coeff[i]==0){
            power--;
                continue;}
        if(power==0){
                if(coeff[i]>0)cout<<"+";
            //else cout<<"-";
               cout<<coeff[i];
        continue;
            }
        if(!m){
            m=true;
            cout<<coeff[i]<<"X^"<<power;
        }
        else{
            if(coeff[i]>0)cout<<"+";
            //else cout<<"-";
            cout<<coeff[i]<<"X^"<<power;
        }
        power--;
     }
     cout<<"=0"<<endl;

 }
 double f(double x){
     double val=0;
     double power=coeff.size()-1;
     for(int i=0;i<coeff.size();i++){
          val+= coeff[i]* pow(x,power);
          power--;
     }
     return val;

 }
int main(){
int degree;
cout<<"enter the degree: ";
cin>>degree;
cout<<"enter the coefficient :";
for(int i=0;i<=degree;i++){
        int x;
      cin>>x;
     coeff.push_back(x);
}
cout<<"equation is: ";
print(coeff);
cout<<endl;
double xmax=0,e=0.001;
for(int i=0;i<coeff.size();i++){
    double temp=coeff[i]/coeff[0];
    xmax=max(xmax,temp);
}
xmax++;
double c=-xmax;
while(c<=xmax){
    double x0=c;
    double x1=c+0.45;
    double fx0=f(x0),fx1=f(x1);
    if(fx0*fx1<0){
        cout<<"search interval is: ["<<x0<<","<<x1<<"]"<<endl;
        int it=0;
        do{
            double x2=x1-f(x1)*((x1-x0)/(fx1-fx0));
            it++;
            double fx2=f(x2);
            if(abs(x1-x2)<e && fx2<e){
                cout<<"root is ="<<x2<<endl<<"iteration is = "<<it<<endl<<endl;
                break;
            }
            x0=x1;
            x1=x2;
            fx0=fx1;
            fx1=fx2;

        }
        while(1);
    }
    c=c+0.45;
}


return 0;
}

```

---

### Newton Raphson

**File:** `NonLinear/newton_raphson.cpp`

```
#include<bits/stdc++.h>
using namespace std;
int degree;
vector<double>coeff;
void print(){
    int d=degree;
    bool m=true;
    for(int i=0;i<degree+1;i++){
            if(m){
            cout<<coeff[i]<<"X^"<<d;
            m=false;
            d--;
            continue;

        }
        if(coeff[i]==0){
                    d--;
                continue;
            }

        if(coeff[i]>0)cout<<"+";
        if(i!=degree)cout<<coeff[i]<<"X^"<<d;
        else cout<<coeff[i];
        d--;
    }
    cout<<"=0"<<endl;
}
 double f(double x){

   double val=0;
     double deg=coeff.size()-1;
     for(int i=0;i<coeff.size();i++){
        val+= coeff[i]* pow(x,deg);
        deg--;
     }
     return val;
 }
 double df(double x){
    double val=0;
     double deg=coeff.size()-1;
     for(int i=0;i<coeff.size()-1;i++){
        val+=coeff[i] *deg* pow(x,deg-1);
        deg--;
     }
     return val;

 }
int main(){
    cout<<"enter degree:";
     cin>>degree;
     cout<<"enter the coefficient:";
     coeff.resize(degree+1);
     for(int i=0;i<degree+1;i++){
        cin>>coeff[i];
     }
     cout<<"equation is: ";
     print();
      double xmax=sqrt((coeff[1]/coeff[0])*(coeff[1]/coeff[0])-2*(coeff[2]/coeff[0]));
      double c=-xmax;
      double num=0,e=0.00001;
      while(c<=xmax){
        double x0=c,x1=c+0.65;
        double fx0=f(x0),fx1=f(x1);
        if(fx0*fx1<0){
                num++;
            cout<<"search interval is ["<<x0<<","<<x1<<"]"<<endl;
            int it=0;
            do{
                double df1=df(x1);
                double x2=x1-(fx1/df1);
              double fx2=f(x2);
              it++;
              if(abs(x2-x1)<e|| abs(fx2-fx1)<e){
                    cout<<"the root "<< num<<" is: "<<x2<<endl;
                 cout<<"iteration is:"<<it<<endl<<endl;
                break;
              }
              x1=x2;
              fx1=fx2;

            }while(1);
        }
        c=c+0.65;
      }

return 0;
}

```

---

## Interpolation

### Newton Forward

**File:** `Interpolation/newton_forward.cpp`  
**Note:** Requires equally spaced `x`.

```
#include <bits/stdc++.h>
using namespace std;

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n; cin >> n;
    vector<double> x(n), y(n);
    for(int i=0;i<n;i++) cin >> x[i] >> y[i];
    double xp; cin >> xp;

    vector<vector<double>> diff(n, vector<double>(n,0.0));
    for(int i=0;i<n;i++) diff[i]=y[i];
    for(int j=1;j<n;j++)
        for(int i=0;i<n-j;i++)
            diff[i][j]=diff[i+1][j-1]-diff[i][j-1];

    double h = x-x;[3]
    double u = (xp - x)/h;

    double yp = diff;
    double uterm = 1.0, fact = 1.0;
    for(int j=1;j<n;j++){
        uterm *= (u-(j-1));
        fact *= j;
        yp += (uterm/fact)*diff[j];
    }

    cout.setf(std::ios::fixed); cout<<setprecision(10);
    cout << "y("<<xp<<") ~ " << yp << "\n";
    return 0;
}
```

---

### Newton Backward

**File:** `Interpolation/newton_backward.cpp`  
**Note:** Requires equally spaced `x`.

```
#include <bits/stdc++.h>
using namespace std;

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n; cin >> n;
    vector<double> x(n), y(n);
    for(int i=0;i<n;i++) cin >> x[i] >> y[i];
    double xp; cin >> xp;

    vector<vector<double>> diff(n, vector<double>(n,0.0));
    for(int i=0;i<n;i++) diff[i]=y[i];
    for(int j=1;j<n;j++)
        for(int i=n-1;i>=j;i--)
            diff[i][j]=diff[i][j-1]-diff[i-1][j-1];

    double h = x-x;[3]
    double u = (xp - x[n-1])/h;

    double yp = diff[n-1];
    double uterm = 1.0, fact = 1.0;
    for(int j=1;j<n;j++){
        uterm *= (u+(j-1));
        fact *= j;
        yp += (uterm/fact)*diff[n-1][j];
    }

    cout.setf(std::ios::fixed); cout<<setprecision(10);
    cout << "y("<<xp<<") ~ " << yp << "\n";
    return 0;
}
```

---

### Divided Difference

**File:** `Interpolation/divided_difference.cpp`  
**Note:** Works for non-equally spaced points.

```
#include <bits/stdc++.h>
using namespace std;

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n; cin >> n;
    vector<double> x(n);
    vector<vector<double>> dd(n, vector<double>(n,0.0));

    for(int i=0;i<n;i++){
        cin >> x[i] >> dd[i];
    }
    double xp; cin >> xp;

    for(int j=1;j<n;j++){
        for(int i=0;i<n-j;i++){
            dd[i][j] = (dd[i+1][j-1]-dd[i][j-1])/(x[i+j]-x[i]);
        }
    }

    double yp = dd;
    double term = 1.0;
    for(int j=1;j<n;j++){
        term *= (xp - x[j-1]);
        yp += term * dd[j];
    }

    cout.setf(std::ios::fixed); cout<<setprecision(10);
    cout << "y("<<xp<<") ~ " << yp << "\n";
    return 0;
}
```

---

## Numerical Integration

> Edit `f(x)` inside each file.



---

## Solution of Numerical Integrations

### Simpson's One-Third Rule

#### Simpson's One-Third Rule Theory
#### Objective
To approximate the definite integral of a function using parabolic interpolation.

#### Data Requirement
A polynomial equation of `n` degree:

  - aₙxⁿ + aₙ₋₁xⁿ⁻¹ + … + a₂x² + a₁x + a₀ = 0

Integration limits: `[a, b]` and number of subintervals `n` (must be even)

#### Core Idea
Fit a parabola → to approximate the curve   
Integrate the parabola → to approximate the area

Divides the interval into an even number of subintervals and fits parabolas through groups of three consecutive points. Once a parabola is fitted through three points, that parabola can be integrated exactly using basic calculus.
**Formula:**
```
∫ₐᵇ f(x)dx ≈ (h/3)[f(x₀) + 4Σf(xᵢ) + 2Σf(xᵢ) + f(xₙ)]
              (i=odd)    (i=even)
```
where h = (b-a)/n

#### Simpson's One-Third Rule Code
```cpp
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
```

```

#### Simpson's One-Third Rule Input
```
3
1 2 0 1
0 1 6
```

#### Simpson's One-Third Rule Output
```
Equation: 1*x^3 + 2*x^2 + 0*x^1 + 1*x^0
Simpson 1/3 Rule Result = 1.91667
```
#### [Back to Contents](#table-of-contents)
---

### Simpson's Three-Eighths Rule 

#### Simpson's Three-Eighths Rule Theory
#### Objective
To approximate the definite integral of a function using cubic interpolation.

#### Data Requirement
A polynomial equation of `n` degree:

  - aₙxⁿ + aₙ₋₁xⁿ⁻¹ + … + a₂x² + a₁x + a₀ = 0

Integration limits: `[a, b]` and number of subintervals `n` (must be multiple of 3)

#### Core Idea
Divides the interval into subintervals (multiple of 3) and fits cubic polynomials through groups of four consecutive points. Provides higher accuracy for functions with higher-order derivatives.  
**Formula:**
```
∫ₐᵇ f(x)dx ≈ (3h/8)[f(x₀) + 3Σf(xᵢ) + 2Σf(xᵢ) + f(xₙ)]
               (i≠3k)      (i=3k)
```
where h = (b-a)/n

#### Simpson's Three-Eighths Rule Code
```cpp
#include <bits/stdc++.h>
using namespace std;

// Function to integrate
float func(float x) {
    float s = sin(x);
    return pow(s, 5) + 4 * pow(s, 4) + 1;
}

int main() {
    float a, b;
    int n;
    cout << "Enter the values of a, b and n (n must be divisible by 3): ";
    cin >> a >> b >> n;
    float h = (b - a) / n;
    float sum = func(a) + func(b);

    for (int i = 1; i < n; i++) {
        float x = a + i * h;
        if (i % 3 == 0) {
            sum += 2 * func(x);
        } else {
            sum += 3 * func(x);
        }
    }

    float result = (3 * h / 8) * sum;
    cout << fixed << setprecision(6);
    cout << "Integral result = " << result << endl;

    return 0;
}
```

#### Simpson's Three-Eighths Rule Input
```
4
1 0 0 1 1
0 3 6
```

#### Simpson's Three-Eighths Rule Output
```
Equation: 1*x^4 + 0*x^3 + 0*x^2 + 1*x^1 + 1*x^0
Simpson 3/8 Rule Result = 56.1562
```



---

## Numerical Differentiation

**File:** `Differentiation/numerical_diff.cpp`  
Computes approximate \(f'(x)\) using forward/backward/central difference.

```
#include <bits/stdc++.h>
using namespace std;

double f(double x){
    return sin(x); // example
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    double x,h;
    cin >> x >> h;

    double forward  = (f(x+h)-f(x))/h;
    double backward = (f(x)-f(x-h))/h;
    double central  = (f(x+h)-f(x-h))/(2*h);

    cout.setf(std::ios::fixed); cout<<setprecision(10);
    cout << "Forward  f'(x) ~ " << forward  << "\n";
    cout << "Backward f'(x) ~ " << backward << "\n";
    cout << "Central  f'(x) ~ " << central  << "\n";
    return 0;
}
```

---

## Least Squares

### Linear Regression

**Model:** \(y = a + bx\)  
**File:** `LeastSquares/linear_regression.cpp`

```
#include <bits/stdc++.h>
using namespace std;
// Linear Regression: y = a + b*x
int main()
{
    int n;
    cout << "Enter the number of data: ";
    cin >> n;

    vector<double> x(n), y(n);

    cout << "Enter the values of x: ";
    for (int i = 0; i < n; i++)
    {
        cin >> x[i];
    }

    cout << "Enter the values of y: ";
    for (int i = 0; i < n; i++)
    {
        cin >> y[i];
    }

    double sx = 0, sy = 0, sxy = 0, sx_2 = 0;

    for (int i = 0; i < n; i++)
    {
        sx += x[i];
        sy += y[i];
        sxy += x[i] * y[i];
        sx_2 += x[i] * x[i];
    }
    double b = (n * sxy - sx * sy) / (n * sx_2 - sx * sx);
    double a = (sy - b * sx) / n;
    cout << fixed << setprecision(4);
    cout << "The linear equation is: y = " << a << " + " << b << "*x" << endl;
    double x_val = 5;
    cout << "The value of y at x = " << x_val << " is: " << (a + b * x_val) << endl;
    return 0;
}

```

---

### Polynomial Regression

**Model:** \(y = a_0 + a_1 x + \dots + a_m x^m\)  
**File:** `LeastSquares/polynomial_regression.cpp`

```
#include <bits/stdc++.h>
using namespace std;

void gaus(vector<vector<double>> &matrix, vector<double> &res)
{
    int n = matrix.size();
    for (int i = 0; i < n; i++)
    {
        int piv = i;
        for (int j = i + 1; j < n; j++)
        {
            if (fabs(matrix[j][i]) > fabs(matrix[piv][i]))
            {
                piv = j;
            }
        }
        swap(matrix[i], matrix[piv]);

        for (int j = i + 1; j < n; j++)
        {
            double r = matrix[j][i] / matrix[i][i];
            for (int k = i; k <= n; k++)
            {
                matrix[j][k] -= r * matrix[i][k];
            }
        }
    }
    res.assign(n, 0);
    for (int i = n - 1; i >= 0; i--)
    {
        res[i] = matrix[i][n];
        for (int j = i + 1; j < n; j++)
        {
            res[i] -= matrix[i][j] * res[j];
        }
        res[i] /= matrix[i][i];
    }
}

int main()
{
    int d, n;
    cout << "Enter the degree: ";
    cin >> d;
    cout << "Enter the number of data: ";
    cin >> n;

    vector<double> x(n), y(n);
    cout << "Enter the value of x: ";
    for (int i = 0; i < n; i++)
        cin >> x[i];

    cout << "Enter the value of y: ";
    for (int i = 0; i < n; i++)
        cin >> y[i];

    vector<vector<double>> matrix(d + 1, vector<double>(d + 2, 0));
    for (int i = 0; i <= d; i++)
    {
        for (int j = 0; j <= d; j++)
        {
            for (int k = 0; k < n; k++)
            {
                matrix[i][j] += pow(x[k], i + j);
            }
        }
        for (int k = 0; k < n; k++)
        {
            matrix[i][d + 1] += y[k] * pow(x[k], i);
        }
    }

    vector<double> res;
    gaus(matrix, res);

    cout << "Polynomial coefficients :";
    for (auto x : res)
        cout << x << " ";
    cout << endl;

    return 0;
}

```

---

### Transcendental Regression

Supports common transforms:
- Exponential: \(y = a e^{bx}\)  →  \(\ln y = \ln a + bx\)
- Power: \(y = a x^b\)          →  \(\ln y = \ln a + b\ln x\)

**File:** `LeastSquares/transcendental_regression.cpp`

```
#include <bits/stdc++.h>
using namespace std;

// fitting a linear line
// y = a + be^(x/4)
// y = a + bx' = a+bf(x)
// x' = f(x) = e^(x/4)

pair<double, double> linear_reg(vector<double> vx, vector<double> vy)
{
    double n = vx.size();
    double x = 0, xy = 0, y = 0, x2 = 0;
    for (int i = 0; i < n; i++)
    {
        x += exp(vx[i] / 4);
        y += vy[i];
        xy += exp(vx[i] / 4) * vy[i];
        x2 += exp(vx[i] / 4) * exp(vx[i] / 4);
    }
    double b = (n * xy - x * y) / (n * x2 - x * x);
    double a = (y - b * x) / n;
    return {b, a};
}

int main()
{
    int n;
    cout << "Enter the no of points: ";
    cin >> n;
    vector<double> vx(n), vy(n);
    cout << "Enter values of x: " << endl;
    for (int i = 0; i < n; i++)
    {
        cin >> vx[i];
        cin >> vy[i];
    }

    auto [b, a] = linear_reg(vx, vy);
    cout << "y = " << a << " + " << b << "e^(x/4)" << endl;
    cout << "Enter any value of x: ";
    double x;
    cin >> x;
    cout << "y(" << x << ") = " << a + b * exp(x / 4) << endl;
}

```

---

## Runge Kutta

**Goal:** Solve \(y' = f(x,y)\) using RK4.  
**File:** `ODE/runge_kutta_rk4.cpp`

```
#include <bits/stdc++.h>
using namespace std;

double f(double x, double y){
    // Example: y' = x + y
    return x + y;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    double x0, y0, h;
    int n;
    cin >> x0 >> y0 >> h >> n;

    double x=x0, y=y0;
    cout.setf(std::ios::fixed); cout<<setprecision(10);
    cout << "x y\n";
    cout << x << " " << y << "\n";

    for(int i=1;i<=n;i++){
        double k1 = h * f(x, y);
        double k2 = h * f(x + h/2.0, y + k1/2.0);
        double k3 = h * f(x + h/2.0, y + k2/2.0);
        double k4 = h * f(x + h,     y + k3);
        y += (k1 + 2*k2 + 2*k3 + k4)/6.0;
        x += h;
        cout << x << " " << y << "\n";
    }
    return 0;
}
```

---
```













