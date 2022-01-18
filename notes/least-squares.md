<!-- @def reeval = true -->
@def hascode = true
\newcommand{\LS}{_\mathrm{\tiny LS}}

# Least Squares

\blurb{Linear least-squares, also known as ordinary linear regression, is a  basic optimization problem used for data fitting. It's also a fundamental building block for algorithms meant for more general problems.}

### Required Reading

- [Sections 3.1--3.2](https://doi.org/10.1137/1.9781611973655.ch3) of Beck. 

## Example: Multilinear regression

We begin with a simple example and use linear regression to build a model that relates the fuel efficiency of a car with its weight, engine displacement, and engine horsepower.  Load the [`mtcars`](https://vincentarelbundock.github.io/Rdatasets/doc/datasets/mtcars.html) dataset provided in the [`Rdatasets.jl`](https://github.com/JuliaStats/RDatasets.jl) package:

```julia:./code/ls1-1
using RDatasets
mtcars = RDatasets.dataset("datasets", "mtcars")
```

```julia:./code/ls-size
# hideall 
print("$(size(mtcars,1))")
```

Here are the first five rows of the dataset, which contains \textoutput{./code/ls-size} rows in total:
```julia:./code/ls1-1
first(mtcars, 5)
```
\show{./code/ls1-1}

The model that we wish to build aims to find coefficients $x_1,\ldots,x_4\in\R$ such that for each particular car model $i=1,2,\ldots,\textoutput{./code/ls-size}$,
\begin{equation*}
\text{MPG}_i \approx x_1 + x_2\cdot\text{Disp}_i + x_3\cdot\text{HP}_i + x_4\cdot\text{WT}_i.
\end{equation*}

Of course, we shouldn't expect an exact match for all car models $i$, and instead the coefficients $x_1,\ldots,x_4$ computed by ordinary linear regression minimizes the sum of least-squares deviations:

\begin{equation} \label{mpg-ls}
\min_{x_1,\ldots,x_4}\ \sum_{i=1}^{\textoutput{./code/ls-size}} (x_1 + x_2\cdot\text{Disp}_i + x_3\cdot\text{HP}_i + x_4\cdot\text{WT}_i - \text{MPG}_i)^2
\end{equation}

Because this is such an important problem in general, there exists excellent software packages that make formulating and solving these problems easy. Here's how to do it using the generalized linear models package [GLM.jl](https://juliastats.org/GLM.jl/stable/): 

```julia:./code/ls-glm
using GLM
ols = lm(@formula(MPG ~ Disp + HP + WT), mtcars)
```
\show{./code/ls-glm}

The column marked `Coef.` shows the values for the 4 computed coefficients $x_1,\ldots,x_4$; the remaining columns show the results for various statistical tests that indicate the reliability of each coefficient
It appears that the car's weight has the biggest impact on fuel economy. This plot shows the projection of the data and the regression line onto the plane spanned by `WT` and `MPG`:

```julia:./code/ls-fit
using StatsPlots; pyplot()
x₁, x₄ = coef(ols)[[1,4]]
@df mtcars scatter(:WT, :MPG, xlabel="weight (1000s lbs)", ylabel="MPG", leg=false)
@df mtcars plot!(:WT, x₁ .+ x₄*:WT)
savefig(joinpath(@OUTPUT,"ls-fit")) # hide
```
\fig{ls-fit}

### Matrix formulation

The objective \eqref{mpg-ls} can be written using matrix notation, as follows: define the matrix and vector
\begin{equation*}
A = \begin{bmatrix}
  1 & \text{WT}_1 & \text{Disp}_1 & \text{HP}_1
\\\vdots & \vdots & \vdots & \vdots
\\1 & \text{WT}_m & \text{Disp}_m & \text{HP}_m
\end{bmatrix},
\quad
b = \begin{bmatrix} \text{MPG}_1 \\ \vdots \\ \text{MPG}_m \end{bmatrix},
\quad
x = \begin{bmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{bmatrix}, 
\end{equation*}
where $m=32$. Thus, the least-squares soluition $x^*$ minimizes the 2-norm squared $\|r\|^2$ of the residual $r = b - Ax$.

\notation{The 2-norm of a vector $z\in\R^n$ is the function
\begin{equation*}
\|z\|_2 = \sqrt{z_1^2+\cdots+z_n^2}.
\end{equation*}
This norm has the property that

\begin{equation*}
  z\T z = \begin{bmatrix} z_1\\\vdots\\ z_n\end{bmatrix}\T
	  \begin{bmatrix} z_1\\\vdots\\ z_n\end{bmatrix}
	= \sum_{i=1}^n z_i^2 = \|z\|_2^2.
\end{equation*}

We often leave out the subscript and simply write $\|z\|$ for the 2-norm, unless warned otherwise.
}

<!-- ## Linear systems

Consider the problem of solving the linear system of equations
$$ Ax = b, $$
where $A \in \R^{m\times n}$ and $b\in\R^{m}$. Then we say this system is

 1. overdetermined if $m>n$,
 2. underdetermined if $m< n$, or
 3. square if $m = n$.

A linear system can have exactly one solution, many solutions, or no solutions:

\fig{./leastsquares_datafitting.png}

But always, a linear system $Ax=b$ has a solution if and only if $b \in \range(A)$. -->

## General least-squares formulation 

The general linear-least squares problem is formulated as
\begin{equation} \label{eq:least-squares-problems}
\min_{x\in\R^n}\ \tfrac12\|r\|^2, \quad r = b - Ax,
\end{equation}
where $A$ is an $m$-by-$n$ matrix, $b$ is an $m$-vector, $x$ is an $n$-vector, and $r$ is an $m$-vector of residuals. Typically, there are many more observations than $m$ than variables $n$, and so $m\gg n$. But there are many important application that require least-squares problems where the opposite is true -- these problems are considered **under determined** because there are necessarily infinitely many least-squares solutions $x\LS.$

The following fundamental result describes several equivalent conditions that any least-squares solution must satisfy.

\theorem{Least-squares optimality}{The vector $x\LS$ solves the least-squares problem \eqref{eq:least-squares-problems} if and only if it satisfies the following equivalent conditions hold:

1. $A^T r = 0$, where $r=b-Ax$,
2. $A\T A x = A^T b$,
3. $y:=Ax$ is the unique orthogonal projection of $b$ onto $\range(A)$.

Moreover, $x\LS$ is the unique least-squares solution if and only if $A$ has full rank.
}

This result relies on a fundamental property that for all matrices $A$, regardless of their shape, the subspaces

\begin{align*}
   \range(A) &= \{y \mid y = A x \text{ for some } x \}
 \\\text{null}(A\T\,) &= \{z \mid A\T z = 0 \}
\end{align*}

are orthogonal complements, i.e.,

\begin{equation*} \label{least_squares_FTLA}
  \range(A) \oplus \Null(A^\intercal) = \R^m.
\end{equation*}

The equivalence of the three properties is then straightforward, and it's sufficient to show that the third condition is equivalent to the optimality of the least-squares problem. Orthogonally decompose the residual $r = b-y$ as

\begin{equation*}
  b - y = z_n + z_r, \quad z_r\in\range(A), \quad z_n\in\Null(A\T).
\end{equation*}

If $y$ was _not_ the orthogonal projection of $b$, then there must exist some component of $b-y$ in the nullspace of $A\T$. Thus, $z_n\ne0$ and

\begin{align*}
  \|r\|^2 &= \|b-y\|^2
\\        &= \|z_n+z_r\|^2
\\        &= \|z_n\|^2+\|z_r\|^2
\\        &> \|z_r\|^2=\|b-(y+z_r)\|^2,
\end{align*}

which means we found a new point $(y+z_r)\in\range(A)$ with a smaller residual, which contradicts the optimality of $x\LS$, and also the uniqueness of the projection $y$.

## Example: Multilinear regression (continued)

For the [multilinear regression example](#example-multilinear-regression), we can obtain the least-squares solution by solving the normal equations as follows:

```julia:./code/ls1
using LinearAlgebra # gives `norm` function
A = [ones(size(mtcars,1)) mtcars[:,:Disp]  mtcars[:,:HP] mtcars[:,:WT]]
b = mtcars[:,:MPG]
x = A'A \ A'b
@show(norm(x-coef(ols)))
```
\show{./code/ls1}
The coefficients computed using the normal equations are _almost_ the same as those computed using by GLM, which likely solves the normal equations using a different method.

The least-squares solution should also verify the first condition of the above theorem:
```julia:./code/ls1
r = b - A*x
@show(norm(A'r))
```
\show{./code/ls1}

## The backslash operator 

In this last example, we solved for the least-squares solution by explicitly solving the normal equations $A\T Ax=A\T b$ via the Julia command `x = A'A \ A'b`. As we'll see in a later lecture, this isn't the best way to solve the least-squares problem. Without getting into the details just yet, it's enough for now to know that Julia allows for the short-hand `x = A \ b`: Julia recognizes that this is an overdetermined system, and solves for the least-squares solution using a mathematically equivalent approach. Let's verify:

```julia:./code/ls1
x1 = A \ b
@show norm(x1 - x, Inf)
```
\show{./code/ls1}

Close enough!

## Example: Trendline 

This example uses least-squares to illustrate the increasing concentrations of [nitrous oxide (N₂O)](https://en.wikipedia.org/wiki/Nitrous_oxide), which is a greenhouse gas. We'll use a dataset from the [Global Monitoring Laboratory](https://gml.noaa.gov):

```julia:./code/ls-n2o
fname = "mlo_N2O_All.dat"
if !isfile(fname)
	download("https://gml.noaa.gov/aftp/data/hats/n2o/insituGCs/CATS/hourly/mlo_N2O_All.dat", fname)
end
```

This dataset provides hourly measurements of N₂O levels in parts per billion (ppb) from 1999 to 2021 from gas sampled at the [Manua Loa Observatory](https://www.google.com/maps/place/Mauna+Loa+Observatory/@19.5363584,-155.5786506,17z/data=!3m1!4b1!4m5!3m4!1s0x7953ef53fcc844c9:0x125658bfa768626b!8m2!3d19.5363852!4d-155.5764517) in Hawaii. The following code removes rows with missing concentration levels `ppb`, and creates a new table with average monthly measurements:

```julia:./code/ls-n2o-1
using CSV, DataFrames, DataFramesMeta, Statistics
df = CSV.read(fname, DataFrame, comment="#", normalizenames=true, delim=" ", ignorerepeated=true)
dfg = @chain df begin
    @select(:year = :N2OcatsMLOyr, :month = :N2OcatsMLOmon, :ppb = :N2OcatsMLOm)
    @subset((!isnan).(:ppb))
    @by([:year, :month], :avgppb = mean(:ppb))
end 
describe(dfg, :min, :max, cols=[:year, :avgppb])
```
\show{./code/ls-n2o-1}

A trend line is defined by the univariate function
\begin{equation*}
  p(t) = \alpha + \beta t,
\end{equation*}
where $\alpha$ is the intercept, $\beta$ is the slope, and $t$ is time. Here, we'll just index the months by the integers $t=1,2,\ldots,T$. We determine the parameters $\alpha$ and $\beta$ by solving the least-squares problem
\begin{equation*}
 \min_{\alpha,\ \beta}\ \sum_{i=1}^{T}(p(t_i) - \text{ppb}_i)^2,
\end{equation*}
which corresponds to the general least-squares problem with
\begin{equation*}
  A = \begin{bmatrix} 1  & t_1 \\ 1 & t_2 \\ \vdots & \vdots \\ 1 & T \end{bmatrix},
  \quad
  b = \begin{bmatrix} \text{ppb}_1 \\ \text{ppb}_2 \\ \vdots \\ \text{ppb}_T \end{bmatrix},
  \textt{and}
  x = \begin{bmatrix} \alpha \\ \beta \end{bmatrix}.
\end{equation*}
In Julia,
```julia:./code/ls-n2o-3
T = nrow(dfg)
t = collect(1:T)
A = [ones(T)  t]
b = dfg[!, :avgppb]
```
The least-squares solution is then
```julia:./code/ls-n2o-4
α, β = A \ b
```
and we can plot the data together with the trend line:
```julia:./code/ls-n2o-2
@df dfg scatter(:avgppb) 
p(t) = α + β*t
plot!(p, lw=3)
plot!(leg=false, xlab="Month index from 1999 to 2021", ylab="Monthly average N₂O ppb")
savefig(joinpath(@OUTPUT,"ls-n2o-2")) # hide
```
\fig{ls-n2o-2}

But perhaps most interestingly, we can now "detrend" the data to understand how the concentrations vary around the overall 22-year trend:
```julia:./code/ls-n2o-3
@df dfg scatter(:avgppb - p.(t))
plot!(leg=false, xlab="Month index from 1999 to 2021", ylab="ppb deviation from mean") # hide
savefig(joinpath(@OUTPUT,"ls-n2o-3")) # hide
```
\fig{ls-n2o-3}

