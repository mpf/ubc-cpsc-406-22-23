---
marp: true
headingDivider: 2
paginate: true
math: true
footer: '[CPSC 406](https://friedlander.io/ubc-cpsc-406)'
---

# Computational Optimization

CPSC 406, Department of Computer Science

Professor Michael P. Friedlander

## Course goals and emphasis

- recognize and formulate the main optimization problem classes
- understand how to apply standard algorithms for each class
- recognize that an algorithm has succeeded or failed
- hands-on experience with mathematical software

## Role of optimization

- fitting a statistical model to data (machine learning)
- logistics, economics, finance, risk management
- theory of games and competition
- theory of computer science and algorithms
- geometry and  analysis

## Learning models

a **model** is a simplified abstraction of process

- model parameters $x=(x_1,x_2,…,x_n)\in\mathcal{P}$
- prediction $m(x)\in\mathcal{F}$
  
least-error principle: the **optimal parameters** $x^*$ minimizes the **distance** between the model and the observation

## Mathematical Optimization

- Objective function: $f:ℝ^n\toℝ$
- feasible set (eg, "constraints"): $\mathcal{C}⊆ℝ^n$
- decision variables: $x=(x_1,x_2,\ldots,x_n)\in\mathcal{C}$

## Abstract problem

- find $x\in\mathcal{C}$ such that $f(x)$ is minimal, eg,
$$
p^* = \min_{x\in\mathcal{C}}\ f(x)
$$
- optimal solution set
$$\mathcal{S}:=\set{x\in\mathcal{C}\mid p^*=f(x)}$$

## Gradients wanted

**Assume** the objective $f$ is differentiable on the **nonempty interior** of  the feasible set $\mathcal{C}\subseteq \mathbb{R}^n$

$$
\nabla f(x) = \begin{pmatrix}
\frac{\partial f(x)}{\partial x_{1}} \\
\vdots\\
\frac{\partial f(x)}{\partial x_{n}}
\end{pmatrix}
$$

- measures the objective's sensitivity to feasible perturbations
- usually sufficient to devise tractable and implementable algorithms

## Example: linear models

$$
\begin{align*}
  \textrm{features}    \qquad A &= [a_1,a_2,\cdots,a_n],\ a_i\inℝ^m
\\\textrm{observations} \qquad b &= (b_1,b_2,\ldots,b_m),\ b\inℝ^m
\\ \textrm{linear model}      \qquad  b &= m(x) + r, \qquad m(x):= Ax
\end{align*}
$$

- least-squares approximation: $x^*=(A^T A)^{-1}A^T b$ minimizes
  $$f(x):=\|r(x)\|^2= \textstyle\sum_i[r(x)]_i^2$$
- least absolute-sum approximation: no closed-form solution for minimizing
  $$f(x):=\|r(x)\|_1 = |r_1(x)| + \cdots + |r_m(x)|$$

## Example: Scheduling

Minimize number nurses needed to meet weekly staffing demands

### Constraints

- each nurse works 5 straight days with 2 days off
- $d_j$ nurses required on nights $j=1,\ldots,7$

### First attempt

- $y_j$ nurses work on night $j$
- minimize $\sum_j y_j$ subject to $y_j\ge d_j$ with $j\in1:7$
- doesn't respect days off constraints

## Scheduling: second attempt

let $x_j$ be number of nurses **starting** their 5-day shift on day $j$:

$$
\begin{array}{ll}
  \min & x_{1} + x_{2} + x_{3} + x_{4} + x_{5} + x_{6} + x_{7} \\
  \text{st}  & \begin{aligned}[t]
           x_{1}\phantom{+x_2+x_3 }\ \ +x_{4}+x_{5}+x_{6}+x_{7}&\ge d_{1}
                 \\x_{1}+x_{2}\phantom{+x_{3}+x_{4}}\ \ +x_{5}+x_{6}+x_{7}&\ge d_{2}
                 \\x_{1}+x_{2}+x_{3}\phantom{+x_{4}+x_{5}}\ \ +x_{6}+x_{7}&\ge d_{3}
                 \\x_{1}+x_{2}+x_{3}+x_{4}\phantom{+x_{5}+x_{6}}\ \ +x_{7}&\ge d_{4}
                 \\x_{1}+x_{2}+x_{3}+x_{4}+x_{5}\phantom{+x_{6}+x_{7}}\ \ &\ge d_{5}
                 \\\phantom{x_{1}+}x_{2}+x_{3}+x_{4}+x_{5}+x_{6}\phantom{+x_{7}}\ \ &\ge d_{6}
                 \\\phantom{x_{1}+x_{2}+}x_{3}+x_{4}+x_{5}+x_{6}+x_{7}&\ge d_{7}
                 \\x_1,\ldots,x_7&\ge0
               \end{aligned}
\end{array}
$$

- note the constraint structure. This is almost always true of practical LPs
- we may want to restrict $x_{j}$ to be integer. That's a much harder problem!

## Coursework and evaluation

- 8 homework assignments (30%)
  - programming and mathematical deriviations
  - typeset submissions, correctness, and writing quality graded
- midterm exam (30%): ✏️️ and 🗞️, short mathematical problems
- final exam (40%): multiple choice

## Programming in Julia

- common programming languages in optimization: [Python](https://www.python.org/), [Matlab](https://www.mathworks.com/), [R](https://www.r-project.org/); C++
- we'll use [Julia](https://julialang.org/) only
  - Matlab-like syntax, but **free** and **fast**
  - lots of [tutorials](https://julialang.org/learning/) available
- [Pluto](https://plutojl.org/) ❤️ or [Jupyter](https://jupyter.org/) or ️notebooks highly recommended for assignments

## Assignments

- work alone or collaborate in pairs 👯
  - hand in your own assignment and list collaborators
- 3 late days allowed (no permission required)
  - **but no more than** 2 late days for a particular assignment
- no solutions posted online
  - come to office hours to see solutions
- submit solutions to Crowdmark

## Homework 1

- no late days, no collaborators
- hints at required background
- due next week

## Resources

- see the [course home page](https://friedlander.io/ubc-cpsc-406) for schedule
- announcements and discussions on [Piazza](https://piazza.com/ubc.ca/winterterm22022/cpsc406)
- TA and instructor office hours start week 2
  