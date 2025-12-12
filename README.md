# **Bisection & False Position Methods**

**Source:** *Lab-2: Bisection and False Position Methods* (CSE 2208 — Numerical Methods Laboratory). 

---

## Overview

* **Course / Lab:** CSE 2208 — Numerical Methods Laboratory.
* **Instructor:** Lamisa Bintee Mizan Deya (Dept. of CSE, KUET).
* **Topic:** Bracketing root-finding methods — **Bisection** and **False Position (Regula Falsi)**. 

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
* Lecturer contact: [deyalamisa49@gmail.com](mailto:deyalamisa49@gmail.com). 

---


* Turn this into a one-page printable note (PDF).
* Expand the worked example with the actual iteration table (I’ll extract the slide values and format them). Which would you prefer?
# NumProject
