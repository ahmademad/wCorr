On the consistent estimation of correlations with sample weights

# Introduction
survey data is weighted to account for varrying inclusion probabilities of individuals and is routinely accounted for in staistics generated with survey data. An exception is the Spearman, polyserial, and polychoric correlations forwhich there is no agreed upon method of generating consistent estimates.

In this paper we present consistent estimators for each correlation for a weighted sample. The main insight is that the nusiance parameters (the districtuion for the polyserial and polychoric, and the ranks themselves) can be weighted to find consistent esimators.

#Defintion of consistent esimation

The correlation   estimators considered in this paper are of the form
$$r_{x,y}(\bm{x},\bm{y}, \bm{w})$$
with $\bm{x} \in R^k$, $\bm{y} \in R^n$ for the Spearman; $\bm{x} \in \{1, 2, ..., m\}$ for some $m<\inf$, $\bm{y} \in R^n$ for the Polyserial; and $\bm{x} \in \{1, 2, ..., m\}$ for some $m<\inf$ $\bm{y} \in \{1, 2, ..., n\}$ for some $n<\inf$ for the polychoric. here $r$ is estimating the opulation value
$$r_{x,y}(\bm{x},\bm{y}, \bm{1})$$
where $k$ is the total number of units in the population and so every unit is self representing, with a weight of one.

Such an estimator is said to be consistent when
$$r_{x,y}(\bm{x},\bm{y}, \bm{w}) \plim \rho$$
as $k \to \inf$ while $\bar{w}=\frac{1}{k} \sum_{i=1}^k w_i$ remains constant. That is, the number of units in the sample and the population grow in proportion to eachother.

#Consistent Estimation of the Spearman Correlations

Estimating a Spearman type correlation requries two steps: first ranks are generated for the two variables, second the Spearman correlation is calculated. When weights are availale, it is not difficult to imagine that the second step should use the formula for the wieghted Pearson correlation.
$$formula for the wieghted Pearson correlation.$$
But this adjustment alone will not result in consistent estimators. If, in addition, the ranks are estimated with weights, so that, when there are no ties,
$$rank(x_i) =\sum_{j=1}^k {\rm 1\kern -2.85 pxI}(x_j \lte x_i) w_j$$
where $I_{ij}$ is an indocator function that is 1 when $x_j \lte x_i$. Notice that because $x_i \lte x_i$, $w_i$ is counted.

When there are ties in at $x_i$, then the following formula is used
$$rank(x_i) =\sum_{j=1}^k {\rm 1\kern -2.85 pxI}(x_j < $x_i) w_j + \frac{n_{tie}+1}{2n_{tie}} \sum_{j=1}^k {\rm 1\kern -2.85 pxI}(x_j = $x_i) w_j$$
where $n_{tie}$ is the number of tied units 
$$n_{tie} = \sum_{j=1}^k {\rm 1\kern -2.85 pxI}(x_j = $x_i)$$
this function has the advantage of returning the typical value of the average rank among the tied units when all of the weights are one while allowing for varying weights.

the estimators is consistent because proof.

#Consistent Estimation of the Polyserial and Polychoric Correlations

The polyserial and polychoric correlations are closely related in that they both accept a ordered covariate for $\bm{x}$ (polyserial) or $\bm{x}$ and $\bm{y}$. The polyserial assumes that there is a latent variable $\bm{x}'$ that is bivariate normally distributed with $\bm{y}$ with correlation coefficient $\rho$. For the polychoric both $\bm{x}$ and $\bm{y} are ordered sets and there exists latent variables $\bm{x}'$ and $\bm{y}'$ that have a bivariate normal distribution with correlation coefficient $\rho$. Both then find $r$ as the MLE.

The variable(s) that are latent are assumed to have cutoffs ($\bm{\zeta}$ for $\bm{x}$ and $\bm{\xi}$ for $\bm{y}$) that define bins so that, for both the polyserial and polychoric,
$$x=i {\rm \ iff} \zeta_{i-1} < y < \zeta_i$$
and also for the polychoric,
$$y=i {\rm \ iff} \xi_{i-1} < y < \xi_i$$
here $\zeta_0 = \xi_0= -\inf$ and $\zeta_m = \xi_n=\inf$, and the $m-1$ finite $\bm{\zeta}$ and $n-1$ finite $\bm{\xi}$ values are nusiance parameters for estimation.

A consistent and ready method of estimating these parameters is 

In the litterature, two methods are popular for estimating $\rho$.

We focus on the follwing method that has better performance in our simulations and is faster.

 and assume that there is a latent variable $\bm{x}'$ (and $\bm{y}'$ in the case of the polychoric) that is jointly bivariately normal

 assume a latent space where one (polyserial) or both (polychoric)

There are two general methods of estimating the polyserial and polychoric correlations


