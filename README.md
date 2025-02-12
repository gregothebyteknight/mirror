## One sphere
Дана сфера фиксированного радиуса R, на основании её делается серия дисекций (например, 10000)
**Целью** является нахождение распределение радиусов (r) серии дисекций
### Ход мысли
Пусть дана сфера с радиусом R и произвольная дисекция с неизвестным радиусом r. Определим расстояние от центра сферы до дисекции как $d=\sqrt{R^2-r²}$ из чего следует, что радиус дисекции равен $r=\sqrt{R² -d²}$
Определимся, что **расстояние d** распределено равномерно и варьирует от -R до R - $d(f)=\frac{1}{R}$ - из соображений симметрии
Отсюда следует, что **распределение r**: $$f(r|R)=f(d)\bullet \frac{\partial d}{\partial r}=\frac{1}{R}\bullet \frac{r}{\sqrt{R^2-r^2}}$$Это же и будет нашим аналитическим решением
Для проверки гипотезы проведём симуляцию проведения 10000 дисекций посредством инициирования расстояний d из равномерного распределения и подсчёта радиусов, дальнейшей визуализации в виде гистограммы
Кроме того, для определения статистического закона  воспользуемся бета распределением (`beta.fit`)
Для сравнения двух выбранных аппроксимаций используем метрику AIC, с кодом для вычисления MLE:
```python
mle_beta = np.sum(beta.logpdf(r, a = a, b = b, loc = loc, scale = scale)) # on the whole r dataset
mle_analytical = np.sum(np.log(analytical_pdf(r, R, eps))) # on the whole r dataset 
```
#### Analytical PDF
Suppose we have function $f(r|R)=\frac{1}{R}\bullet \frac{r}{\sqrt{R^2-r^2}}$ and using the Bayes formula $f(R|r)=\frac{f(r|R)f(R)}{f(r)}$ it can derived that
$$f(r)=\int_{R=r}^{\infty}f(r(|R)f(R)$$
On avera
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
