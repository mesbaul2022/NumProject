# Numerical Lab Project (C++)

Numerical Methods Laboratory implementations in **C++**: solution of linear equations, non-linear root finding, interpolation, numerical integration & differentiation, least squares regression, and Runge–Kutta (ODE).

---

## Index
- [Solution of Linear Equations](#solution-of-linear-equations)
  - [Gauss Elimination](#gauss-elimination)
  - [Gauss Jordan](#gauss-jordan)
  - [LU Decomposition](#lu-decomposition)
  - [Matrix Inversion](#matrix-inversion)
- [Solution of Non-linear Equations](#solution-of-non-linear-equations)
  - [Bisection](#bisection)
  - [False Position](#false-position)
  - [Secant](#secant)
  - [Newton Raphson](#newton-raphson)
- [Interpolation](#interpolation)
  - [Newton Forward](#newton-forward)
  - [Newton Backward](#newton-backward)
  - [Divided Difference](#divided-difference)
- [Numerical Integration](#numerical-integration)
  - [Simpson 1/3](#simpson-13)
  - [Simpson 3/8](#simpson-38)
- [Numerical Differentiation](#numerical-differentiation)
- [Least Squares](#least-squares)
  - [Linear Regression](#linear-regression)
  - [Polynomial Regression](#polynomial-regression)
  - [Transcendental Regression](#transcendental-regression)
- [Runge Kutta](#runge-kutta)

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

### Simpson 1/3

**File:** `Integration/simpson_1_3.cpp`  
**Note:** `n` must be even.

```
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

---

### Simpson 3/8

**File:** `Integration/simpson_3_8.cpp`  
**Note:** `n` must be a multiple of 3.

```
#include <bits/stdc++.h>
using namespace std;

double f(double x){
    return 1.0/(1.0+x*x); // example
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    double a,b; int n;
    cin >> a >> b >> n;
    if(n%3!=0){
        cout << "n must be multiple of 3.\n";
        return 0;
    }

    double h = (b-a)/n;
    double sum = f(a) + f(b);

    for(int i=1;i<n;i++){
        double x = a + i*h;
        sum += (i%3==0 ? 2.0 : 3.0) * f(x);
    }

    double ans = (3.0*h/8.0)*sum;

    cout.setf(std::ios::fixed); cout<<setprecision(10);
    cout << "Integral ~ " << ans << "\n";
    return 0;
}
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

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n; cin >> n;
    vector<double> x(n), y(n);
    for(int i=0;i<n;i++) cin >> x[i] >> y[i];

    double sx=0, sy=0, sxx=0, sxy=0;
    for(int i=0;i<n;i++){
        sx += x[i]; sy += y[i];
        sxx += x[i]*x[i];
        sxy += x[i]*y[i];
    }

    double den = n*sxx - sx*sx;
    if (fabs(den) < 1e-15){
        cout << "Cannot fit (denominator near zero).\n";
        return 0;
    }
    double b = (n*sxy - sx*sy)/den;
    double a = (sy - b*sx)/n;

    cout.setf(std::ios::fixed); cout<<setprecision(10);
    cout << "a = " << a << "\n";
    cout << "b = " << b << "\n";
    cout << "Fitted: y = a + b*x\n";
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

static const double EPS = 1e-12;

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, m;
    cin >> n >> m;              // n points, degree m
    vector<double> x(n), y(n);
    for(int i=0;i<n;i++) cin >> x[i] >> y[i];

    int dim = m+1;
    vector<vector<double>> A(dim, vector<double>(dim,0.0));
    vector<double> B(dim,0.0);

    // Build normal equations
    vector<double> sx(2*m+1,0.0);
    for(int k=0;k<=2*m;k++){
        for(int i=0;i<n;i++) sx[k] += pow(x[i], k);
    }
    for(int i=0;i<dim;i++){
        for(int j=0;j<dim;j++) A[i][j] = sx[i+j];
        for(int t=0;t<n;t++) B[i] += y[t]*pow(x[t], i);
    }

    // Solve A * coeff = B using Gauss elimination
    vector<vector<double>> aug(dim, vector<double>(dim+1));
    for(int i=0;i<dim;i++){
        for(int j=0;j<dim;j++) aug[i][j]=A[i][j];
        aug[i][dim]=B[i];
    }

    for(int col=0; col<dim; col++){
        int pivot=col;
        for(int r=col;r<dim;r++)
            if (fabs(aug[r][col]) > fabs(aug[pivot][col])) pivot=r;
        if (fabs(aug[pivot][col]) < EPS){
            cout << "Cannot fit (singular normal matrix).\n";
            return 0;
        }
        swap(aug[pivot], aug[col]);
        for(int r=col+1;r<dim;r++){
            double factor = aug[r][col]/aug[col][col];
            for(int k=col;k<=dim;k++) aug[r][k]-=factor*aug[col][k];
        }
    }

    vector<double> c(dim);
    for(int i=dim-1;i>=0;i--){
        double sum=aug[i][dim];
        for(int j=i+1;j<dim;j++) sum-=aug[i][j]*c[j];
        c[i]=sum/aug[i][i];
    }

    cout.setf(std::ios::fixed); cout<<setprecision(10);
    for(int i=0;i<dim;i++) cout << "a" << i << " = " << c[i] << "\n";
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

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n; cin >> n;
    int mode; 
    // 1 = exponential (y = a e^(b x))
    // 2 = power       (y = a x^b)
    cin >> mode;

    vector<double> X(n), Y(n);
    for(int i=0;i<n;i++) cin >> X[i] >> Y[i];

    vector<double> u(n), v(n);
    for(int i=0;i<n;i++){
        if (Y[i] <= 0){
            cout << "All y must be > 0 for log transform.\n";
            return 0;
        }
        if (mode==1){
            u[i] = X[i];
            v[i] = log(Y[i]);
        } else if (mode==2){
            if (X[i] <= 0){
                cout << "All x must be > 0 for power model.\n";
                return 0;
            }
            u[i] = log(X[i]);
            v[i] = log(Y[i]);
        } else {
            cout << "Invalid mode.\n";
            return 0;
        }
    }

    double su=0, sv=0, suu=0, suv=0;
    for(int i=0;i<n;i++){
        su += u[i]; sv += v[i];
        suu += u[i]*u[i];
        suv += u[i]*v[i];
    }

    double den = n*suu - su*su;
    if (fabs(den) < 1e-15){
        cout << "Cannot fit (denominator near zero).\n";
        return 0;
    }

    double b = (n*suv - su*sv)/den;
    double A = (sv - b*su)/n; // A = ln(a)
    double a = exp(A);

    cout.setf(std::ios::fixed); cout<<setprecision(10);
    cout << "a = " << a << "\n";
    cout << "b = " << b << "\n";
    if (mode==1) cout << "Fitted: y = a * exp(b*x)\n";
    if (mode==2) cout << "Fitted: y = a * x^b\n";
    return 0;
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

























# **Bisection & False Position Methods**


---

* **Course / Lab:** CSE 2208 — Numerical Methods Laboratory.
* **Topic:** **Bisection** and **False Position (Regula Falsi)**. 

---

## Key concepts

* **Equation types:** linear vs non-linear; root finding focuses on non-linear equations (f(x)=0).
* **Bracketing methods:** start with two guesses (a) and (b) such that (f(a)) and (f(b)) have opposite signs (they *bracket* a root). Methods reduce the bracket until desired accuracy is achieved. 

---

## Bisection Method — Essentials

* **Idea:** Repeatedly bisect interval ([a,b]) and pick subinterval where sign change occurs.
* **Step-by-step:**

  1. Ensure (f(a)\cdot f(b) < 0).
  2. Compute midpoint (c = \tfrac{a+b}{2}).
  3. If (f(c)=0) (or |f(c)| small), stop. Else choose new bracket: if (f(a)f(c)<0) set (b=c), else set (a=c).
  4. Repeat until (|b-a|) or (|f(c)|) meets tolerance.
* **Properties:** Guaranteed convergence (linear), robust, error halves each iteration. Slow compared with open methods. Figures and an example solution table are provided in the lab slides. 

---

## False Position (Regula Falsi) — Essentials

* **Idea:** Use a secant line through ((a,f(a))) and ((b,f(b))); the x-intercept replaces one bracket endpoint.
* **Formula:** (c = b - f(b)\dfrac{b-a}{f(b)-f(a)}).
* **Step-by-step:**

  1. Start with (a,b) s.t. (f(a)f(b)<0).
  2. Compute (c) using the secant intercept formula.
  3. Replace the endpoint (a) or (b) that has the same sign as (f(c)) with (c).
  4. Repeat until tolerance satisfied.
* **Properties:** Usually faster than bisection (uses function values), but can suffer slow convergence if one endpoint remains nearly fixed. Lab includes illustration and an example table (noting a case of only 3 iterations).

---

## Worked example(s) (as in the slides)

* **Test function used in problems:**
  [
  f(x) = 3x - \cos x - 1
  ]
  This function is used for both bisection and false position examples; the slides contain iteration tables and graphs for the solutions. 

---

## Quick comparison (at a glance)

* **Reliability:** Bisection ✅ (always converges) — False Position ✅ (converges but may stagnate).
* **Speed:** False Position usually faster than Bisection when it behaves well.
* **Simplicity:** Both simple; Bisection uses midpoint, False Position uses a secant intercept.
* **When to use:** Bisection for guaranteed, robust root brackets; False Position to accelerate when function behavior is favorable. 

---

## Minimal pseudocode

**Bisection**

```
Given a, b with f(a)*f(b) < 0, tolerance tol
while (b - a)/2 > tol:
    c = (a + b)/2
    if f(c) == 0: return c
    if f(a)*f(c) < 0: b = c
    else: a = c
return (a + b)/2
```

**False Position**

```
Given a, b with f(a)*f(b) < 0, tolerance tol
while |f(c)| > tol:
    c = b - f(b)*(b - a)/(f(b) - f(a))
    if f(a)*f(c) < 0: b = c
    else: a = c
return c
```

---

## Practical notes & tips

* Check continuity and that the function actually changes sign on the initial bracket.
* Monitor both interval width (|b-a|) and function value (|f(c)|) for stopping.
* If False Position stalls (one endpoint doesn’t change), consider modified regula-falsi variants or switch to bisection/secant hybrid. 

---

## References (from the lab)

* YouTube Playlist 1 & 2 (numerical methods videos).
* **Book:** *Numerical Methods* by E. Balagurusamy.

---


* Turn this into a one-page printable note (PDF).
* Expand the worked example with the actual iteration table (I’ll extract the slide values and format them). Which would you prefer?
# NumProject



















# Lab-3: Newton–Raphson & Secant — **CSE 2208: Numerical Methods Laboratory**

---

## 1) Overview — *Open methods (interpolation / iterative)*

* **Definition:** Open methods start with one or two initial guesses (they do **not** require the root to be bracketed) and iteratively refine those guesses using interpolation formulas until convergence.
* **Two methods covered:** **Newton–Raphson** and **Secant**. 

---

## 2) Newton–Raphson method — key facts

* **Iteration formula**
  [
  x_{n+1}=x_n - \frac{f(x_n)}{f'(x_n)}
  ]
* **Convergence:** quadratic (fast) if (f'(x^\star)\neq 0) and initial guess is sufficiently close. Requires computation of derivative (f').
* **Geometric idea:** approximate (f) by its tangent at (x_n), step to the tangent root (see slide figure). 

---

### Worked example from the slides (Problem 1)

**Problem statement (slide):**
Find the root *in the vicinity of* (x=0) of
[
f(x)=3x-\cos x - 1.
]
(Problem & slides: Lab-3). 

* **Derivative**
  [
  f'(x)=3+\sin x.
  ]

* **Newton–Raphson example** (chosen initial guess (x_0=0.5), shown to illustrate the iterations):

|  n |     (x_n)    |    (f(x_n))   |   (f'(x_n))  |   (x_{n+1})  |
| -: | :----------: | :-----------: | :----------: | :----------: |
|  1 | 0.5000000000 | -0.3775825619 | 3.4794255386 | 0.6085186499 |
|  2 | 0.6085186499 |  0.0050602142 | 3.5716526463 | 0.6071018788 |
|  3 | 0.6071018788 | 8.2373687e-07 | 3.5704896184 | 0.6071016481 |
|  4 | 0.6071016481 | 2.1760371e-14 | 3.5704894289 | 0.6071016481 |

* **Converged root (NR):** (x^\star \approx \mathbf{0.6071016481031226}).

*(Iterations above are a compact numerical run of the NR formula to show rapid quadratic convergence.)*

---

## 3) Secant method — key facts

* **Iteration formula**
  [
  x_{n+1}=x_n - f(x_n),\frac{x_n-x_{n-1}}{f(x_n)-f(x_{n-1})}
  ]
  (uses two last iterates; no derivative required).
* **Convergence:** superlinear (≈1.618 order), typically slower than Newton–Raphson but avoids evaluating (f').
* **Geometric idea:** replace the tangent by the secant line through two recent points. See slide figure and Problem 2. 

---

### Worked example from the slides (Problem 2)

**Problem statement (slide):** same function
[
f(x)=3x-\cos x - 1,
]
solved by the **Secant** method. (Slides show the problem & a result table; initial guesses were used in the lecture.) 

* **Secant example** (illustrative initial guesses (x_0=0.0,; x_1=1.0)):

|  n |   (x_{n-1})  |     (x_n)    |  (f(x_{n-1})) |    (f(x_n))    |   (x_{n+1})  |
| -: | :----------: | :----------: | :-----------: | :------------: | :----------: |
|  1 | 0.0000000000 | 1.0000000000 | -2.0000000000 |  1.4596976941  | 0.5780851903 |
|  2 | 1.0000000000 | 0.5780851903 |  1.4596976941 |  -0.1032549064 | 0.6059585719 |
|  3 | 0.5780851903 | 0.6059585719 | -0.1032549064 |  -0.0040808049 | 0.6071055027 |
|  4 | 0.6059585719 | 0.6071055027 | -0.0040808049 |  1.3762708e-05 | 0.6071016476 |
|  5 | 0.6071055027 | 0.6071016476 | 1.3762708e-05 | -1.8100770e-09 | 0.6071016481 |

* **Converged root (Secant):** (x^\star \approx \mathbf{0.6071016481031226}) (same root as NR).

*(This demonstrates Secant converging to the same root; secant requires a couple more iterations here but avoids derivative evaluation.)*

---

## 4) Practical notes — when to use which

* **Newton–Raphson**

  * Pros: very fast (quadratic) near the root.
  * Cons: need (f'); may diverge if initial guess poor or (f'(x)) small.
* **Secant**

  * Pros: no derivative required; cheaper per iteration.
  * Cons: slower than NR (superlinear), needs two initial guesses; may fail if iterates give equal function values.

---

## 5) References (from the provided slides)

* YouTube Playlists (lecture support) — links in the slides. 
* Book: *Numerical Methods* by E. Balagurusamy. 

---


















