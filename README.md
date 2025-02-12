## One sphere
A sphere of fixed radius R is given, and a series of dissections (e.g. 10000) is made based on it
**Objective** is to find the distribution of radii (r) of a series of dissections
### Thought Process 
Let a sphere with radius R and an arbitrary dissection with unknown radius r be given. Let us define the distance from the center of the sphere to the dissection as $d=\sqrt{R^2-r²}$ from which it follows that the radius of the dissection is $r=\sqrt{R² -d²}$

Let us define that **distance d** is uniformly distributed and varies from -R to R - $d(f)=\frac{1}{R}$ - from 0 to R for reasons of symmetry 

It follows that **distribution r**: $$f(r|R)=f(d)\bullet \frac{\partial d}{\partial r}=\frac{1}{R}\bullet \frac{r}{\sqrt{R^2-r^2}}$$. This will be our analytical solution

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
**Plotting R inference boxplots**: `est_box_const()` from *`viz.py`* < `big_r_est()` from *`const_case.py`* < `sim()` from *`const_case.py

## Set of spheres
For this iteration I chose a set of spheres (~1000) with R distributed by a gamma law. For this computational experiment I simulate one dissection with 1000 small radii. 
**Aim** is to infer the distribution of small radii
### Pipeline
For each of R from gamma law I first compute distance and then r. From a set of small radii and plot a histogram which will be a picture of my simulation
As a pattern descriptor I use two metrics: iterative analytical function and fully analytical function
#### Iterative Analytical Function

#### Fully Analytical Function

## R Inference
### Single R Case
This task implies construction of the R estimator based on r distribution. Here, we consider two estimators: derived from **mean of r distribution** and **max of r**
Also we need to check our estimator based on **two criteria**: convergence and absence of bias

#### Mean Estimator
$$\mathbb{E}(f(r))=\mathbb{E}(\frac{1}{R}\frac{r}{\sqrt{R^2 -r^2}})=\frac{1}{R}\int_0^R\frac{r^2}{\sqrt{R^2-r^2}}=\frac{R\pi}{4}$$ Which mean that we can try to use $R=\mathbb{E}(r)\bullet\frac{4}{\pi}$ as our estimator

#### Max Estimator
$$\bar{R}_n=max(r)$$

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
