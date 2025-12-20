Gaussian elimination is commonly described as a step-by-step procedure that solves linear systems by applying row operations to transform the augmented matrix into (upper) row‑echelon form before back‑substitution.[1]
Below is a complete, ready-to-paste `README.md` structure (with all your requested topics) that you can edit to match your exact file names and languages.

```markdown
# Numerical Methods Laboratory

A collection of Numerical Methods lab implementations covering **linear** systems, **non-linear** equations, **interpolation**, **numerical integration/differentiation**, **least squares regression**, and **Runge–Kutta** methods.

---

## Table of Contents

- [Project Overview](#project-overview)
- [Features](#features)
- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
  - [How to Run](#how-to-run)
  - [Input/Output Convention](#inputoutput-convention)
- [Implemented Methods](#implemented-methods)
  - [Linear Methods](#linear-methods)
    - [Gauss Elimination](#gauss-elimination)
    - [Gauss-Jordan](#gauss-jordan)
    - [LU Decomposition](#lu-decomposition)
    - [Matrix Inversion](#matrix-inversion)
  - [Non-linear Methods](#non-linear-methods)
    - [Bisection Method](#bisection-method)
    - [False Position (Regula Falsi)](#false-position-regula-falsi)
    - [Secant Method](#secant-method)
    - [Newton–Raphson Method](#newtonraphson-method)
  - [Interpolation](#interpolation)
    - [Newton Forward Interpolation](#newton-forward-interpolation)
    - [Newton Backward Interpolation](#newton-backward-interpolation)
    - [Newton Divided Difference](#newton-divided-difference)
  - [Numerical Integration](#numerical-integration)
    - [Simpson’s 1/3 Rule](#simpsons-13-rule)
    - [Simpson’s 3/8 Rule](#simpsons-38-rule)
  - [Numerical Differentiation](#numerical-differentiation)
  - [Least Squares Regression](#least-squares-regression)
    - [Linear Regression](#linear-regression)
    - [Polynomial Regression](#polynomial-regression)
    - [Transcendental Regression](#transcendental-regression)
  - [Runge–Kutta (ODE)](#rungekutta-ode)
- [Project Structure](#project-structure)
- [Adding a New Method](#adding-a-new-method)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgements](#acknowledgements)

---

## Project Overview

This repository is designed for Numerical Methods laboratory work and practice.
Each method includes a short description, expected inputs/outputs, and code implementation.

---

## Features

- Clear separation by topic (Linear / Non-linear / Interpolation / Integration / etc.)
- Standard tolerance-based stopping criteria for iterative methods
- Example input format included for quick testing
- Easy-to-extend structure for adding more methods

---

## Getting Started

### Prerequisites

Update this section depending on your language:

- **C/C++**: GCC / G++ (or CodeBlocks)
- **Python**: Python 3.x
- **MATLAB/Octave** (optional): If you used MATLAB/Octave scripts

### How to Run

> Replace file names below with your actual file paths.

#### C/C++ (example)
```
# Compile
g++ path/to/file.cpp -o run

# Run
./run
```

#### Python (example)
```
python path/to/file.py
```

### Input/Output Convention

Use a consistent format across programs (recommended):

- **Inputs** (typical):
  - Matrix size `n`
  - Coefficient matrix `A`
  - RHS vector `b`
  - Initial guesses / bracket `[a, b]`
  - Tolerance `eps`
  - Maximum iteration `maxIter`

- **Outputs** (typical):
  - Approximated solution (vector or root)
  - Iteration table (optional but recommended)
  - Error value per iteration (optional)

---

## Implemented Methods

## Linear Methods

### Gauss Elimination

**Goal:** Solve a system of linear equations \(Ax=b\) using forward elimination + back substitution.

**Inputs**
- `n` (order of the system)
- `A` (coefficient matrix)
- `b` (constant vector)
- (Optional) pivoting on/off

**Outputs**
- Solution vector `x`
- (Optional) transformed upper-triangular matrix

**Algorithm (high-level)**
- Build the augmented matrix \([A|b]\)
- Forward eliminate below the diagonal
- Back substitute to obtain `x`

**Code**
- File: `./Linear/GaussElimination/...`
- Link: _Add your GitHub file link here_

**Sample**
- Input: _Put a small example here_
- Output: _Expected output here_

---

### Gauss-Jordan

**Goal:** Solve \(Ax=b\) by transforming \([A|b]\) to reduced row echelon form (RREF).

**Inputs**
- `n`, `A`, `b`

**Outputs**
- Solution vector `x`
- (Optional) RREF matrix

**Algorithm (high-level)**
- Make pivot 1
- Make all other entries in pivot column 0
- Read solution directly

**Code**
- File: `./Linear/GaussJordan/...`

---

### LU Decomposition

**Goal:** Factorize \(A = LU\) and solve \(Ax=b\) via forward/back substitution.

**Common Variants**
- Doolittle (diagonal of `L` is 1)
- Crout (diagonal of `U` is 1)

**Inputs**
- `n`, `A`, `b`

**Outputs**
- Matrices `L` and `U`
- Solution vector `x`

**Algorithm (high-level)**
- Decompose `A` into `L` and `U`
- Solve `Ly=b` (forward substitution)
- Solve `Ux=y` (back substitution)

**Code**
- File: `./Linear/LU/...`

---

### Matrix Inversion

**Goal:** Compute \(A^{-1}\) (if it exists), typically using Gauss-Jordan on \([A|I]\).

**Inputs**
- `n`, `A`

**Outputs**
- `A_inv` or message: "Singular matrix"

**Algorithm (high-level)**
- Form augmented matrix \([A|I]\)
- Apply Gauss-Jordan
- Extract inverse from the right block

**Code**
- File: `./Linear/MatrixInverse/...`

---

## Non-linear Methods

### Bisection Method

**Goal:** Find a root of \(f(x)=0\) in an interval \([a,b]\) where \(f(a)\cdot f(b) < 0\).

**Inputs**
- Function `f(x)`
- `a`, `b`
- `eps`, `maxIter`

**Outputs**
- Approximated root `x`
- Iteration count (optional)

**Algorithm (high-level)**
- Compute midpoint `c=(a+b)/2`
- Choose subinterval that preserves sign change
- Stop when error < `eps` or `f(c)==0`

**Code**
- File: `./NonLinear/Bisection/...`

---

### False Position (Regula Falsi)

**Goal:** Find a root in \([a,b]\) using a secant line but keeping the bracket.

**Inputs**
- `f(x)`, `a`, `b`, `eps`, `maxIter`

**Outputs**
- Approximated root

**Algorithm (high-level)**
- Compute intersection point using straight line through \((a,f(a))\) and \((b,f(b))\)
- Replace `a` or `b` depending on sign
- Repeat until convergence

**Code**
- File: `./NonLinear/FalsePosition/...`

---

### Secant Method

**Goal:** Find a root using secant updates (no derivative needed).

**Inputs**
- `f(x)`
- Two initial guesses `x0`, `x1`
- `eps`, `maxIter`

**Outputs**
- Approximated root

**Algorithm (high-level)**
- Iteratively compute the next point from the secant line through the last two points

**Code**
- File: `./NonLinear/Secant/...`

---

### NewtonRaphson Method

**Goal:** Find a root using derivative-based iteration.

**Inputs**
- `f(x)`, `f'(x)`
- Initial guess `x0`
- `eps`, `maxIter`

**Outputs**
- Approximated root

**Algorithm (high-level)**
- Update: `x_{k+1} = x_k - f(x_k)/f'(x_k)`
- Stop when error < `eps` or iteration limit reached

**Code**
- File: `./NonLinear/NewtonRaphson/...`

---

## Interpolation

### Newton Forward Interpolation

**Use when:** points are equally spaced and `x` is near the beginning of the table.

**Inputs**
- Data points \((x_i, y_i)\) (equal spacing)
- Target `x`

**Outputs**
- Approximated `y(x)`

**Code**
- File: `./Interpolation/NewtonForward/...`

---

### Newton Backward Interpolation

**Use when:** points are equally spaced and `x` is near the end of the table.

**Inputs**
- Equally spaced \((x_i, y_i)\)
- Target `x`

**Outputs**
- Approximated `y(x)`

**Code**
- File: `./Interpolation/NewtonBackward/...`

---

### Newton Divided Difference

**Use when:** x values are NOT equally spaced (works for general spacing).

**Inputs**
- General \((x_i, y_i)\)
- Target `x`

**Outputs**
- Approximated `y(x)`

**Code**
- File: `./Interpolation/DividedDifference/...`

---

## Numerical Integration

### Simpson's 1/3 Rule

**Goal:** Approximate \(\int_a^b f(x)\,dx\) using parabolic segments.

**Inputs**
- `f(x)`, `a`, `b`
- `n` (must be even)

**Outputs**
- Approximate integral value

**Code**
- File: `./Integration/Simpson_1_3/...`

---

### Simpson's 3/8 Rule

**Goal:** Approximate \(\int_a^b f(x)\,dx\) using cubic interpolation segments.

**Inputs**
- `f(x)`, `a`, `b`
- `n` (must be multiple of 3)

**Outputs**
- Approximate integral value

**Code**
- File: `./Integration/Simpson_3_8/...`

---

## Numerical Differentiation

**Goal:** Approximate derivatives when symbolic differentiation is difficult.

**Typical Methods Included**
- Forward difference
- Backward difference
- Central difference
- Higher-order formulas (if implemented)

**Inputs**
- `f(x)` or tabulated points
- Point `x`
- Step size `h`

**Outputs**
- Approximation of \(f'(x)\) (and optionally \(f''(x)\))

**Code**
- File: `./Differentiation/...`

---

## Least Squares Regression

### Linear Regression

**Model**
- Fit: \(y = a + bx\)

**Inputs**
- Data points \((x_i, y_i)\)

**Outputs**
- Parameters `a`, `b`
- (Optional) predicted values and error metrics

**Code**
- File: `./LeastSquares/Linear/...`

---

### Polynomial Regression

**Model**
- Fit: \(y = a_0 + a_1 x + a_2 x^2 + \dots + a_m x^m\)

**Inputs**
- \((x_i, y_i)\), degree `m`

**Outputs**
- Coefficients \(a_0..a_m\)

**Code**
- File: `./LeastSquares/Polynomial/...`

---

### Transcendental Regression

**Common Models**
- Exponential: \(y = a e^{bx}\)
- Power: \(y = a x^b\)
- Logarithmic: \(y = a + b\ln(x)\)

**Inputs**
- \((x_i, y_i)\)
- Selected model type

**Outputs**
- Model parameters and fitted equation

**Code**
- File: `./LeastSquares/Transcendental/...`

---

## RungeKutta (ODE)

**Goal:** Numerically solve an ODE like \(y' = f(x,y)\) with an initial condition \(y(x_0)=y_0\).

**Common Version**
- 4th order Runge–Kutta (RK4)

**Inputs**
- Function `f(x, y)`
- `x0`, `y0`
- Step size `h`
- Number of steps `n` (or final `x`)

**Outputs**
- Approximated values of `y` over the interval
- (Optional) table of `(x, y)` values

**Code**
- File: `./ODE/RungeKutta/...`

---

## Project Structure

> Edit folder names to match your repository.

```
Numerical-Methods-Laboratory/
│
├─ Linear/
│  ├─ GaussElimination/
│  ├─ GaussJordan/
│  ├─ LU/
│  └─ MatrixInverse/
│
├─ NonLinear/
│  ├─ Bisection/
│  ├─ FalsePosition/
│  ├─ Secant/
│  └─ NewtonRaphson/
│
├─ Interpolation/
│  ├─ NewtonForward/
│  ├─ NewtonBackward/
│  └─ DividedDifference/
│
├─ Integration/
│  ├─ Simpson_1_3/
│  └─ Simpson_3_8/
│
├─ Differentiation/
│
├─ LeastSquares/
│  ├─ Linear/
│  ├─ Polynomial/
│  └─ Transcendental/
│
├─ ODE/
│  └─ RungeKutta/
│
└─ README.md
```

---

## Adding a New Method

1. Create a folder under the correct category.
2. Add code file(s) and (optional) sample input/output.
3. Add a short section under **Implemented Methods** with:
   - Goal
   - Inputs/Outputs
   - Algorithm (high-level)
   - File link/path

---

## Contributing

Contributions are welcome.

- Fork the repository
- Create a feature branch
- Commit changes with clear messages
- Open a pull request

---

## License

Add a license if you want (recommended). Example:

This project is licensed under the MIT License - see the `LICENSE` file for details.

---

## Acknowledgements

- Course: Numerical Methods Laboratory
- Instructor/Department: _Add your info_
- Teammates/Contributors: _Add names (optional)_
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


















