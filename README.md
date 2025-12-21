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
  - [Gauss Jordan Elimination Method](#gauss-jordan-elimination-method)
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
  - [Simpson's One-Third Rule](#simpsons-one-third-rule)
    - [Theory](#simpsons-one-third-rule-theory)
    - [Code](#simpsons-one-third-rule-code)
    - [Input](#simpsons-one-third-rule-input)
    - [Output](#simpsons-one-third-rule-output)
  - [Simpson's Three-Eighths Rule](#simpsons-three-eighths-rule)
    - [Theory](#simpsons-three-eighths-rule-theory)
    - [Code](#simpsons-three-eighths-rule-code)
    - [Input](#simpsons-three-eighths-rule-input)
    - [Output](#simpsons-three-eighths-rule-output)
  


---

## Solution of Linear Equations

### Gauss Elimination

#### Gauss Elimination Theory

#### Method used
**Gauss Elimination Method**

#### Objective
To solve a system of linear algebraic equations by transforming it into an **upper triangular system**, followed by **back substitution**.

#### Data Requirement
A system of `n` linear equations:
```
a₁₁x₁ + a₁₂x₂ + ... + a₁ₙxₙ = b₁
a₂₁x₁ + a₂₂x₂ + ... + a₂ₙxₙ = b₂
...
aₙ₁x₁ + aₙ₂x₂ + ... + aₙₙxₙ = bₙ
```

Matrix form:
```
AX = B
```

#### Notation
- `A = [aᵢⱼ]` : coefficient matrix of order `n × n`
- `X = [x₁, x₂, ..., xₙ]ᵀ` : vector of unknowns
- `B = [b₁, b₂, ..., bₙ]ᵀ` : constant vector

#### Core Idea
The system is simplified by eliminating variables using **elementary row operations** to obtain an upper triangular matrix.

#### Elimination Approach (Formula)
To eliminate element `aᵢⱼ` (where `j < i`):
```
Rᵢ ← Rᵢ − (aᵢⱼ / aⱼⱼ) Rⱼ
```


#### Phases Involved

##### Forward Elimination
Transforms the augmented matrix `[A | B]` into an **upper triangular form**.

##### Back Substitution
Solutions are obtained using:
```
xₙ = bₙ / aₙₙ
xᵢ = (1 / aᵢᵢ) [ bᵢ − Σ (aᵢⱼ xⱼ) ], j = i+1 to n
```
#### Gauss Elimination Code
```cpp

#include <bits/stdc++.h>
using namespace std;

const float EPS = 1e-6;

int main() {
    int t;
    cout << "Enter number of test cases: ";
    cin >> t;

    for (int i = 1; i <= t; i++) {
        int n;
        cout << "Enter number of variables " << i << ": ";
        cin >> n;

        vector<vector<float>> matrix(n, vector<float>(n + 1));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= n; j++) {
                cin >> matrix[i][j];
            }
        }

        for (int i = 0; i < n; i++) {
            int pivotrow = i;
            while (pivotrow < n && abs(matrix[pivotrow][i]) < EPS) {
                pivotrow++;
            }

            if (pivotrow == n) continue;
            swap(matrix[i], matrix[pivotrow]);

            float pivot = matrix[i][i];
            for (int j = i; j <= n; j++) {
                matrix[i][j] /= pivot;
            }

            for (int k = 0; k < n; k++) {
                if (k != i) {
                    float factor = matrix[k][i];
                    for (int j = i; j <= n; j++) {
                        matrix[k][j] -= factor * matrix[i][j];
                    }
                }
            }
        }
        bool nosolution = false;
        bool infinitesolution = false;
        int rank = 0;

        for (int i = 0; i < n; i++) {
            bool allzero = true;
            for (int j = 0; j < n; j++) {
                if (abs(matrix[i][j]) > EPS) {
                    allzero = false;
                    break;
                }
            }

            if (allzero) {
                if (abs(matrix[i][n]) > EPS) {
                    nosolution = true;
                    break;
                } else {
                    infinitesolution = true;
                }
            } else {
                rank++;
            }
        }

        cout << "Case " << i << ":" << endl;
        if (nosolution) {
            cout << "No Solution" << endl;
        } else if (infinitesolution || rank < n) {
            cout << "Infinite Solutions" << endl;
        } else {
            cout << "Unique Solution:" << endl;
            for (int i = 0; i < n; i++) {
                cout<< matrix[i][n] <<" ";
            }
        }
        cout << endl;

    }

    return 0;
}

```

#### Gauss Elimination Input
```
3
3
1 1 1 6
0 2 5 -4
2 5 -1 27
3
1 2 3 5
2 4 6 12
1 1 1 4
3
1 2 3 5
2 4 6 10
1 1 1 4

```

#### Gauss Elimination Output
```
Case 1:
Unique Solution:
5 3 -2

Case 2:
No Solution

Case 3:
Infinite Solutions
```
#### [Back to Contents](#table-of-contents)
---

### Gauss Jordan Elimination Method

#### Gauss Jordan Theory
#### Method used
**Gauss–Jordan Elimination Method**

#### Objective
To solve a system of linear equations by reducing the augmented matrix directly to **Reduced Row Echelon Form (RREF)**.

#### Data Requirement
Augmented matrix form:
```
[A | B]
```

#### Notation
- `aᵢⱼ` : element of coefficient matrix
- `bᵢ`  : element of constant vector

#### Core Idea
Each pivot element is made **unity**, and all other elements in its column are eliminated, producing an identity matrix.

#### Elimination Approach (Formulae)

**Normalization of pivot row:**
```
Rᵢ ← Rᵢ / aᵢᵢ
```

**Elimination of other rows:**
```
Rⱼ ← Rⱼ − aⱼᵢ Rᵢ (j ≠ i)
```
#### Matrix Form Obtained
```
[I | X]
```

where `I` is the identity matrix and `X` contains the solution.

#### Gauss Jordan Code
```cpp
#include <bits/stdc++.h>
using namespace std;

void solve() {
    int n;
    cout<<"Enter the number of variable ";
    cin >> n;
    vector<vector<double>> a(n, vector<double>(n + 1));
    double EPS = 1e-12;
    cout<<"Enter the matrix "<<endl;
    for (int i = 0; i < n; i++)
        for (int j = 0; j <= n; j++) cin >> a[i][j];

    for (int i = 0; i < n; i++) {
        int pivotrow = i;
            while (pivotrow < n && abs(a[pivotrow][i]) < EPS) {
                pivotrow++;
            }

            if (pivotrow == n) continue;
            swap(a[i], a[pivotrow]);

        double t = a[i][i];
        for (int j = 0; j <= n; j++) a[i][j] /= t;

        for (int j = 0; j < n; j++) {
            if (i != j) {
                double r = a[j][i];
                for (int k = 0; k <= n; k++) a[j][k] -= r * a[i][k];
            }
        }
    }

    bool inf = false, none = false;
    for (int i = 0; i < n; i++) {
        bool zero = true;
        for (int j = 0; j < n; j++) if (fabs(a[i][j]) > EPS) zero = false;
        if (zero) {
            if (fabs(a[i][n]) > EPS) none = true;
            else inf = true;
        }
    }

    if (none) cout << "NO solution" << endl;
    else if (inf) cout << "INFINITE Solution" << endl;
    else {
        cout << "Unique solution: ";
        for (int i = 0; i < n; i++) cout << a[i][n] << " ";
        cout << endl;
    }
}

int main() {
    int t;
    cout << "Enter number of test cases: ";
    cin >> t;
    while (t--) {
        solve();
        cout<<endl;
    }
    return 0;
}

```

#### Gauss Jordan Input
```
3
3
1 1 1 6
0 2 5 -4
2 5 -1 27
3
1 2 3 5
2 4 6 12
1 1 1 4
3
1 2 3 5
2 4 6 10
1 1 1 4
```

#### Gauss Jordan Output
```
Unique solution: 5 3 -2

NO solution

INFINITE Solution

```
#### [Back to Contents](#table-of-contents)
---

### LU Decomposition Method

#### LU Decomposition Theory
#### Method used
**LU Decomposition Method**

#### Objective
To solve a system of linear equations by factorizing the coefficient matrix into lower and upper triangular matrices.

#### Data Requirement
A square matrix with non-zero pivots.

#### Core Idea (Formula)
```
A = LU
```

#### Notation
- `L = [lᵢⱼ]` : lower triangular matrix
- `U = [uᵢⱼ]` : upper triangular matrix
- `A` : coefficient matrix
- `X` : solution vector
- `B` : constant vector

#### Solution Process

**Step 1:** Solve
```
LY = B
```
using forward substitution:
```
yᵢ = bᵢ − Σ (lᵢⱼ yⱼ), j = 1 to i−1
```

**Step 2:** Solve
```
UX = Y
```
using back substitution:
```
xᵢ = (1 / uᵢᵢ) [ yᵢ − Σ (uᵢⱼ xⱼ) ], j = i+1 to n
```

#### Accuracy Considerations
- More efficient than repeated Gauss Elimination
- Numerical stability improves with pivoting

#### Applicability
- Ideal for solving multiple systems with the same coefficient matrix


#### LU Decomposition Code
```cpp
#include <bits/stdc++.h>
using namespace std;

int main()
{

    // ---- File Handling ----
    ifstream in("input.txt");
    ofstream out("output.txt");

    if (!in)
    {
        cerr << "Error: input.txt not found\n";
        return 1;
    }

    int casing;
    in >> casing;
    int r = 1;

    while (casing--)
    {
        cout << endl;
        out << "----- Case " << r++ << " -----\n\n";
        int n;
        in >> n;

        vector<vector<double>> a(n + 1, vector<double>(n + 2));

        // Read augmented matrix A|b
        for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= n + 1; j++)
            {
                in >> a[i][j];
            }
        }

        vector<vector<double>> u(n + 1, vector<double>(n + 1, 0));
        vector<vector<double>> l(n + 1, vector<double>(n + 1, 0));

        for (int i = 1; i <= n; i++)
        {
            l[i][i] = 1;
        }

        // ---- LU Decomposition ----
        for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= n; j++)
            {

                if (i <= j)
                {
                    u[i][j] = a[i][j];
                    for (int k = 1; k < i; k++)
                        u[i][j] -= l[i][k] * u[k][j];
                }
                else
                {
                    l[i][j] = a[i][j];
                    for (int k = 1; k < j; k++)
                        l[i][j] -= l[i][k] * u[k][j];

                    if (u[j][j] == 0)
                    {
                        out << "Matrix is singular. Cannot compute LU decomposition.\n";
                        return 0;
                    }

                    l[i][j] /= u[j][j];
                }
            }
        }

        // ---- Print U Matrix ----
        out << "U Matrix:\n";
        for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                out << u[i][j] << " ";
            }
            out << "\n";
        }

        // ---- Print L Matrix ----
        out << "\nL Matrix:\n";
        for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                out << l[i][j] << " ";
            }
            out << "\n";
        }

        // ---- Forward Substitution: Ly = b ----
        vector<double> y(n + 1, 0);

        for (int i = 1; i <= n; i++)
        {
            y[i] = a[i][n + 1];
            for (int k = 1; k < i; k++)
            {
                y[i] -= l[i][k] * y[k];
            }
        }

        // ---- Check Solution Type ----
        bool noSolution = false;
        bool infiniteSolution = false;

        for (int i = 1; i <= n; i++)
        {
            bool allZero = true;

            for (int j = 1; j <= n; j++)
            {
                if (fabs(u[i][j]) > 1e-9)
                {
                    allZero = false;
                    break;
                }
            }

            if (allZero)
            {
                if (fabs(y[i]) > 1e-9)
                {
                    noSolution = true;
                }
                else
                {
                    infiniteSolution = true;
                }
            }
        }

        if (noSolution)
        {
            out << "\nThe system has NO SOLUTION (Inconsistent equations).\n";
            continue;
        }

        if (infiniteSolution)
        {
            out << "\nThe system has INFINITE SOLUTIONS (Dependent equations).\n";
            continue;
        }

        out << "\nThe system has a UNIQUE SOLUTION.\n";

        // ---- Backward Substitution: Ux = y ----
        vector<double> ans(n + 1, 0);

        for (int i = n; i >= 1; i--)
        {
            ans[i] = y[i];
            for (int k = i + 1; k <= n; k++)
            {
                ans[i] -= u[i][k] * ans[k];
            }
            ans[i] /= u[i][i];
        }

        // ---- Print Solution ----
        out << "\nFinal Solution (x values):\n";
        for (int i = 1; i <= n; i++)
        {
            out << "x" << i << " = " << ans[i] << "\n";
        }
    }

    return 0;
}


```

#### LU Decomposition Input
```
4
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3
2
1 1 2
2 2 4
2
1 1 2
2 2 5
5
2 1 -1 3 2 9
1 3 2 -1 1 8
3 2 4 1 -2 20
2 1 3 2 1 17
1 -1 2 3 4 15
```

#### LU Decomposition Output
```
----- Case 1 -----

U Matrix:
2 1 -1 
0 0.5 0.5 
0 0 -1 

L Matrix:
1 0 0 
-1.5 1 0 
-1 4 1 

The system has a UNIQUE SOLUTION.

Final Solution (x values):
x1 = 2
x2 = 3
x3 = -1
----- Case 2 -----

U Matrix:
1 1 
0 0 

L Matrix:
1 0 
2 1 

The system has INFINITE SOLUTIONS (Dependent equations).
----- Case 3 -----

U Matrix:
1 1 
0 0 

L Matrix:
1 0 
2 1 

The system has NO SOLUTION (Inconsistent equations).
----- Case 4 -----

U Matrix:
2 1 -1 3 2 
0 2.5 2.5 -2.5 0 
0 0 5 -3 -5 
0 0 0 1.4 3 
0 0 0 0 1.85714 

L Matrix:
1 0 0 0 0 
0.5 1 0 0 0 
1.5 0.2 1 0 0 
1 0 0.8 1 0 
0.5 -0.6 0.8 1.71429 1 

The system has a UNIQUE SOLUTION.

Final Solution (x values):
x1 = 5.15385
x2 = -1
x3 = 2.26154
x4 = -0.138462
x5 = 1.18462


```
#### [Back to Contents](#table-of-contents)
---

### Matrix Inversion

#### Matrix Inversion Theory

#### Method used
Cofactor/Adjugate-based inverse with recursive determinant.

#### Objective
Solve `AX = B` via `X = A⁻¹B` when `det(A) ≠ 0`; classify singular cases when det(A) = 0.

#### Data Requirement
Square matrix A. If `det(A) = 0`, the system may be inconsistent (no solution) or dependent (infinitely many); classification uses the augmented matrix.

#### Core Idea (Formula)
```
AX = B
```
```
X = A⁻¹ B (if det(A) ≠ 0)
```

#### Notation
- `A⁻¹`: inverse of A
- `adj(A)`: adjugate of A
- `Cᵀ`: transpose of cofactor matrix C

#### Inversion Approach
- Compute det(A) using cofactors.
- If det(A) = 0:
    - If a row of A is all zeros but the corresponding B entry is nonzero → No solution (inconsistent).
    - If rows of A are dependent and B is compatible → Infinite solutions (dependent).
    - Stop.
- Build cofactor matrix C with:
```
Cᵢⱼ = (−1)^(i+j) · det(Mᵢⱼ)
```
- Form adjugate and inverse:
```
adj(A) = Cᵀ
A⁻¹ = adj(A) / det(A)
```

#### Evaluation Process
Multiply A⁻¹ by B to obtain the solution vector X.
```
X = A⁻¹ B
```

#### Accuracy Considerations
- Recursive determinant and cofactors are computationally expensive and numerically sensitive.
- Prefer elimination methods for large systems; this approach is illustrative and aligns with the provided code.

#### Applicability
- Good for educational purposes and small systems where clarity of the inverse construction is desired.


#### Matrix Inversion Code
```cpp
#include <bits/stdc++.h>
using namespace std;

// Cofactor
void getCofactor(const vector<vector<double>>& A,
                 vector<vector<double>>& temp,
                 int p, int q, int n)
{
    int i = 1, j = 1;
    for (int row = 1; row <= n; row++) {
        for (int col = 1; col <= n; col++) {
            if (row != p && col != q) {
                temp[i][j++] = A[row][col];
                if (j == n) {
                    j = 1;
                    i++;
                }
            }
        }
    }
}

//Recursive Determinant
double determinant(const vector<vector<double>>& A, int n)
{
    if (n == 1)
        return A[1][1];

    double det = 0;
    int sign = 1;
    vector<vector<double>> temp(n, vector<double>(n, 0));

    for (int i = 1; i <= n; i++) {
        getCofactor(A, temp, 1, i, n);
        det += sign * A[1][i] * determinant(temp, n - 1);
        sign = -sign;
    }
    return det;
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if (!fin) {
        cout << "Error opening input file.\n";
        return 1;
    }

    int casing;
    fin>> casing;
    int r=1;
    while(casing--){
    fout<<endl;
    fout<<"Case "<<r++<<":\n";
    int n;
    fin >> n;

    vector<vector<double>> aug(n+1, vector<double>(n+2, 0));
    vector<vector<double>> a(n+1, vector<double>(n+1, 0));
    vector<vector<double>> B(n+1, vector<double>(2, 0));
    vector<vector<double>> C(n+1, vector<double>(n+1, 0));
    vector<vector<double>> C1(n+1, vector<double>(n+1, 0));
    vector<vector<double>> res(n+1, vector<double>(2, 0));

    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n + 1; j++)
            fin >> aug[i][j];

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++)
            a[i][j] = aug[i][j];
        B[i][1] = aug[i][n+1];
    }

    double detA = determinant(a, n);

    if (fabs(detA) < 1e-9) {
        // Row reduce to check rank by gauss jordan elimination
        vector<vector<double>> tempAug = aug;
        const double EPS = 1e-9;
        
        for (int i = 1; i <= n; i++) {
            // Find pivot
            int maxRow = i;
            for (int k = i + 1; k <= n; k++) {
                if (fabs(tempAug[k][i]) > fabs(tempAug[maxRow][i]))
                    maxRow = k;
            }
            swap(tempAug[i], tempAug[maxRow]);
            
            if (fabs(tempAug[i][i]) < EPS) continue;
            
            // Eliminate below
            for (int k = i + 1; k <= n; k++) {
                double factor = tempAug[k][i] / tempAug[i][i];
                for (int j = i; j <= n + 1; j++) {
                    tempAug[k][j] -= factor * tempAug[i][j];
                }
            }
        }
        
        // Check for inconsistency:  row with all zeros in A but non-zero in b
        bool noSol = false;
        for (int i = 1; i <= n; i++) {
            bool allZero = true;
            for (int j = 1; j <= n; j++) {
                if (fabs(tempAug[i][j]) > EPS) {
                    allZero = false;
                    break;
                }
            }
            if (allZero && fabs(tempAug[i][n+1]) > EPS) {
                noSol = true;
                break;
            }
        }

        if (noSol)
            fout << "Determinant = 0 → No Solution (Inconsistent System)\n";
        else
            fout << "Determinant = 0 → Infinite Solutions (Dependent System)\n";

        continue;
    }

    fout << "Determinant = " << detA << "\n\n";

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            vector<vector<double>> temp(n, vector<double>(n, 0));
            getCofactor(a, temp, i, j, n);
            C[i][j] = pow(-1, i + j) * determinant(temp, n - 1);
        }
    }

    fout << "Inverse Matrix:\n";
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            C1[i][j] = C[j][i] / detA;
            fout << C1[i][j] << " ";
        }
        fout << "\n";
    }

    for (int i = 1; i <= n; i++) {
        for (int t = 1; t <= n; t++)
            res[i][1] += C1[i][t] * B[t][1];
    }

    fout << "\nSolution Vector:\n";
    for (int i = 1; i <= n; i++)
        fout << "x" << i << " = " << res[i][1] << "\n";
    fout << "\n";
}

    fin.close();
    fout.close();


    return 0;
}

```

#### Matrix Inversion Input
```
4
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3
2
1 1 2
2 2 4
2
1 1 2
2 2 5
5
2 1 -1 3 2 9
1 3 2 -1 1 8
3 2 4 1 -2 20
2 1 3 2 1 17
1 -1 2 3 4 15
```

#### Matrix Inversion Output
```

Case 1:
Determinant = -1

Inverse Matrix:
4 3 -1 
-2 -2 1 
5 4 -1 

Solution Vector:
x1 = 2
x2 = 3
x3 = -1


Case 2:
Determinant = 0 → Infinite Solutions (Dependent System)

Case 3:
Determinant = 0 → No Solution (Inconsistent System)

Case 4:
Determinant = 65

Inverse Matrix:
0.384615 0.384615 1.92308 -3.76923 1.61538 
-0 0 -1 2 -1 
-0.246154 -0.0461538 -0.230769 0.692308 -0.153846 
-0.0461538 -0.446154 -1.23077 2.69231 -1.15385 
0.0615385 0.261538 0.307692 -0.923077 0.538462 

Solution Vector:
x1 = 5.15385
x2 = -1
x3 = 2.26154
x4 = -0.138462
x5 = 1.18462


```
#### [Back to Contents](#table-of-contents)
---

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
#### [Back to Contents](#table-of-contents)
---
## Download PDF

- Download this README as a PDF for offline viewing: [README.pdf](README.pdf)


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

## Solution of Curve Fitting Model

### Least Square Regression Method For Linear Equations Method

#### Least Square Regression Method For Linear Equations Theory
#### Method used
**Least Squares Regression – Linear Equation**

#### Objective
To fit a straight line of the form 
```
y = a + bx
```
that best represents the given experimental data.

#### Data Requirement
A set of `n` observed data points:
```
(x1, y1), (x2, y2), ..., (xn, yn)
```

#### Core Idea
The best-fitting straight line is obtained by minimizing the **sum of squares of vertical deviations** (errors) between the observed data points and the assumed line.

#### Assumed Model
```
y = a + bx
```

#### Least Squares Conditions
To minimize error:
```
∑(y − a − bx)² → minimum
```
This leads to the **normal equations**:
```
∑y = na + b∑x
∑xy = a∑x + b∑x²
```

#### Least Square Regression Method For Linear Equations Code
```cpp
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

#### Least Square Regression Method For Linear Equations Input
```
5
1 2 3 4 5
2 4 5 4 5
```

#### Least Square Regression Method For Linear Equations Output
```
Enter the number of data: 5
Enter the values of x: 1 2 3 4 5
Enter the values of y: 2 4 5 4 5
The linear equation is: y = 2.2000 + 0.6000*x
The value of y at x = 5.0000 is: 5.2000

```
#### [Back to Contents](#table-of-contents)
---

### Least Square Regression Method For Transcendental Equations 

#### Least Square Regression Method For Transcendental Equations Theory
#### Method used
**Least Squares Regression – Transcendental Equation**

#### Objective
To fit nonlinear relationships by transforming them into linear form so that least squares method can be applied.

#### Common Transcendental Forms
1. **Exponential model**
```
y = ae^(bx)
```

2. **Power model**
```
y = ax^b
```

#### Linearization Technique

##### Exponential Equation
Taking natural logarithm:
```
ln y = ln a + bx
```
Let:
```
Y = ln y , A = ln a
```
Then:
```
Y = A + bx
```

##### Power Equation
Taking logarithm on both sides:
```
log y = log a + b log x
```

Let:
```
Y = log y , X = log x , A = log a
```
Then:
```
Y = A + bX
```
#### Evaluation Process
- Transform the given data
- Apply linear least squares regression
- Compute constants
- Convert back to original form

#### Accuracy Considerations
- Transformation may amplify errors
- Requires positive data values
- Fit quality depends on correct model assumption

#### Applicability
- Used in population growth, decay processes, and empirical laws
- Suitable for nonlinear experimental data


#### Least Square Regression Method For Transcendental Equations Code
```cpp
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

#### Least Square Regression Method For Transcendental Equations Input
```
5
0 2
1 3.718
2 6.389
3 10.109
4 15.0
2
```

#### Least Square Regression Method For Transcendental Equations Output
```
Enter the no of points: 5
Enter values of x:
0 2
1 3.718
2 6.389
3 10.109
4 15.0
y = -5.94885 + 7.63686e^(x/4)
Enter any value of x: 2
y(2) = 6.64221

```
#### [Back to Contents](#table-of-contents)
---

### Least Square Regression Method For Polynomial Equations 

#### Least Square Regression Method For Polynomial Equations Theory
#### Method used
**Least Squares Regression – Polynomial Equation**

#### Objective
To fit a polynomial curve when linear regression is insufficient to represent data trends.

#### Assumed Model (Second Order Polynomial)
```
y = a + bx + cx²
```

#### Data Requirement
A set of experimental observations:
```
(x1, y1), (x2, y2), ..., (xn, yn)
```

#### Core Idea
The coefficients of the polynomial are determined by minimizing the sum of squared deviations between observed and computed values.

#### Normal Equations
For a second-degree polynomial:
```
∑y = na + b∑x + c∑x²
∑xy = a∑x + b∑x² + c∑x³
∑x²y = a∑x² + b∑x³ + c∑x⁴
```
#### Evaluation Process
- Compute required summations from data table
- Solve the system of equations
- Substitute coefficients into polynomial

#### Accuracy Considerations
- Higher degree improves fit but may cause overfitting
- Computational complexity increases with degree
- Sensitive to data errors

#### Applicability
- Used when data shows curvature
- Suitable for engineering and experimental modeling

#### Least Square Regression Method For Polynomial Equations Code
```cpp
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

#### Least Square Regression Method For Polynomial Equations Input
```
2
5
1 2 3 4 5
1 4 9 16 25

```

#### Least Square Regression Method For Polynomial Equations Output
```
Enter the degree: 2
Enter the number of data: 5
Enter the value of x: 1 2 3 4 5
Enter the value of y: 1 4 9 16 25
Polynomial coefficients :0 0 1 

```
#### [Back to Contents](#table-of-contents)
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













