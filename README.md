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































***

# CSE 2208 — Numerical Methods Laboratory (Lab 2)[1]

A short lab on **bracketing methods** for finding roots of non-linear equations, focusing on Bisection and False Position methods.[1]

## Table of Contents
- Overview
- Methods Covered
- Lab Problems
- References & Resources

## Overview
This lab is part of **CSE 2208 Numerical Methods Laboratory**.[1]
The slide deck credits Lamisa Bintee Mizan Deya (Lecturer, Dept. of CSE, KUET) and includes an email contact: `deyalamisa49gmail.com`.[1]
It introduces equations (including linear vs non-linear) and positions root-finding as an iterative approach that starts with two initial guesses that bracket a root and then shrinks that bracket.[1]

## Methods Covered
The lab discusses two bracketing-based root-finding techniques: **Bisection Method** and **False Position Method**.[1]
It includes illustrations for both methods (figures shown in the slides) to demonstrate how intervals/brackets update across iterations.[1]

## Lab Problems
- **Problem 1 (Bisection):** Find the root of \( f(x) = 3x\cos(x) - 1 \) using the bisection method.[1]
- **Problem 2 (Bisection continuation):** A second bisection-based solution table is included as “Table 2” in the slides (problem statement/details are not fully visible in the extracted text).[1]
- **Problem 3 (False Position):** Find the root of \( f(x) = 3x\cos(x) - 1 \) using the false position method; the slides note it takes “only 3 iterations this time” (as presented in the deck’s solution table).[1]

## References & Resources
The slide deck points to the following learning resources.[1]
- YouTube Playlist 1: https://youtube.com/playlist?list=PLgH5QX0i9K3oKFrSOo4Kwns1-vTZmKQ7zsiWKta7GwWfg5gAfAE[1]
- YouTube Playlist 2: https://youtube.com/playlist?list=PLU6SqdYcYsfIk1VhXxIYNPFU67ym6gae8siDf0wuN88YiBJFwwI[1]
- Book: *Numerical Methods* by E. Balagurusamy.[1]

--- 



