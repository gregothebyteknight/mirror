## One sphere
A sphere of fixed radius R is given, and a series of dissections (e.g. 10000) is made based on it
**Objective** is to find the distribution of radii (r) of a series of dissections
### Thought Process 
Let a sphere with radius R and an arbitrary dissection with unknown radius r be given. Let us define the distance from the center of the sphere to the dissection as $d=\sqrt{R^2-r²}$ from which it follows that the radius of the dissection is $r=\sqrt{R² -d²}$

Let us define that **distance d** is uniformly distributed and varies from -R to R - $d(f)=\frac{1}{R}$ - from 0 to R for reasons of symmetry 

It follows that **distribution r**: $$f(r|R)=f(d)\bullet \frac{\partial d}{\partial r}=\frac{1}{R}\bullet \frac{r}{\sqrt{R^2-r^2}}$$ This will be our analytical solution

In order to test the hypothesis, we will conduct a simulation of 10000 dissections by initiating the distances d from a uniform distribution and counting the radii, further visualization in the form of a histogram

In addition, we will use the beta distribution to determine the statistical law  (`beta.fit`)

As a second approximation, a separate function (`analytical_pdf`) was written which used the Bayes rule to determine the true marginal distribution of r.

To compare the two selected approximations, we use the AIC metric, with the code to compute MLE:
```python
mle_beta = np.sum(beta.logpdf(r, a = a, b = b, loc = loc, scale = scale)) # on the whole r dataset
mle_analytical = np.sum(np.log(analytical_pdf(r, R, eps))) # on the whole r dataset 
```
#### Analytical PDF
Suppose we have function $f(r|R)=\frac{1}{R}\bullet \frac{r}{\sqrt{R^2-r^2}}$ and using the Bayes formula $f(R|r)=\frac{f(r|R)f(R)}{f(r)}$ it can derived that
$$f(r)=\int_{R=r}^{\infty}f(r(|R)f(R)$$
Approach is presented in the `analytical_pdf` function from *`modules.py`*

#### Visualization Pipelines
**Plotting simulations**: `pdf_const`  from *`viz.py`* < `sim()` and `analytical_pdf` from *`modules.py`*

**Plotting R inference boxplots**: `est_box_const()` from *`viz.py`* < `big_r_est()` from *`const_case.py`* < `sim()` from *`const_case.py`*

## Set of spheres
For this iteration I chose a set of spheres (~1000) with R distributed by a gamma law. For this computational experiment I simulate one dissection with 1000 small radii. 

**Aim** is to infer the distribution of small radii

### Pipeline
For each of R from gamma law I first compute distance and then r. From a set of small radii and plot a histogram which will be a picture of my simulation
As a pattern descriptor I use **two metrics**: iterative analytical function and fully analytical function
#### Iterative Analytical Function
Used as a first attempt to visualize the PDF. Calculation based on the following steps:
1. 

#### Fully Analytical Function

## R Inference
### Single R Case
This task implies construction of the R estimator based on r distribution. Here, we consider two estimators: derived from **mean of r distribution** and **max of r**

Also we need to check our estimator based on **two criteria**: convergence and absence of bias

#### Mean Estimator
$$\mathbb{E}(f(r))=\mathbb{E}(\frac{1}{R}\frac{r}{\sqrt{R^2 -r^2}})=\frac{1}{R}\int_0^R\frac{r^2}{\sqrt{R^2-r^2}}=\frac{R\pi}{4}$$ 
From this formula we can postulate that $R=\mathbb{E}(r)\bullet\frac{4}{\pi}$ and use it as our estimator

#### Max Estimator
$$\bar{R}_n=max(r)$$
Was not used for the simulations

#### Code for Single R Case
Is presented in the `est_box_const` function which for each number of simulations from the list provided runs the r simulation function `n_rep` times. then computes the original R with mean estimator. Calculations the are averaged by the `n_rep` times, statistics of the R inference finally is visualized with boxplots.

Clearly we can see that our mean estimator 

### R Distribution Case
**Conditions**: 1000 follicles with radii $R \sim \gamma(\alpha, \beta)$, one dissection with 1000 radii $r \sim \gamma'(\alpha, \beta)$

**Problem**: as for each of follicles we have only one observation, using simple formula $R=\mathbb{E}(r)\bullet \frac{4}{\pi}$ will give us very noisy estimations as we will not be able to use the true mean for our formula

Here the **two approaches** will be used:
- Deriving parameters from the **equation** $\int_{R=r}^{\infty}\frac{1}{R}\frac{r}{\sqrt{R^2-r^2}}\gamma(\alpha,\beta)dR=\gamma'(\alpha',\beta')$
- Using **assumption**: $r\sim \gamma'(f(\alpha), g(\beta))$
#### Using Assumptions
$$
\left\{\begin{matrix}
R \sim \gamma(\alpha, \beta) \\
r \sim f(R)
\end{matrix}\right.
$$Given above we can assume $r \sim \gamma'(f(\alpha), g(\beta))=\gamma'(\alpha',\beta')$
If ${r_i}\sim \gamma'(\alpha',\beta')$, then we can suggest that $$\left\{\begin{matrix}
\alpha = f^{-1}(\alpha') \\
\beta = g^{-1}(\beta')
\end{matrix}\right.$$**Steps**:
1. Based on 30 repetitions plot the real parameter of gamma (from initialized R gamma) against shifted one (r gamma)
2. Fit the curve with different assumption functions
3. Compute the inverse
4. Implement the inverse into the R inference function 

# Pair-Correlation Function

## Quick Introduction
Pair-correlation function $g(r)$ measures the probability of finding a pair of objects separated by a distance r relative to what would be expected for a completely random (Poisson) distribution
### Derivation
Given that space is isotropic let suggest that probability of observing a point around location x is $$p(x)=\lambda(x)dx$$where $\lambda(x)$ is intensity, or density of the point process (like Poisson process) and $dV$ - small region of volume

Let $p(x,y)$ be a probability of observing 1 point in x and 1 point in y, then $$p(x,y)=\lambda(x)dx\lambda(y)dy\bullet g(x,y)$$where $g(x,y)=\frac{p(x,y)}{\lambda^2dxdy}$ - **pair-correlation function** (in simple case)

### Intuition
Layers of spheres get more diffuse, so for large distances, the probability of finding two spheres with a given separation is essentially constant > a more dense system has more spheres, this it's more likely to find two of them with a given distance

The PC function accounts for these factors by normalizing by the density; this at large values of r it goes to 1, uniform probability

### Calculation in 3D Case
1. Pick a value of `dr`
2. Loop over all values of r
    1. Count all particles that are a distance between `r` and `r+dr` away from the particle you're considering
    2. Divide your total count by N > average local density
    3. Divide this number by $4\pi r^2dr$ (volume of the spherical shell)
    4. Divide this by the particle number density $\rho$ - ensures that $g(r)=1$

For 2D volume correction will be just $2\pi r dr$
