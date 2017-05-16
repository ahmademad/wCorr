---
title: "wCorr Proofs"
author: "Paul Bailey"
date: '`r Sys.Date()`'
header-includes:
  - \usepackage{bm}
  - \usepackage{amsmath, amssymb}
output: pdf_document
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{wCorr Formulas}
  \usepackage[utf8]{inputenc}
---

# Proof of consistency of Horvitz-Thompson estimator of a mean.

An Horvitz-Thompson (HT) estimator of a sum takes the form
\begin{align} \label{sum1}
\hat{Y} &= \sum_{i=1}^n \frac{1}{\pi_i} y_i \ ,
\end{align}
where there are $n$ sampled units from a population of $N$ units, each unit has a value $y \in \Re$, each unit is sampled with probability $\pi_i$, and  $\hat{Y}$ is the estimated of $Y$ in a population. Probabilisticaly, the formula is
\begin{align} \label{sum2}
\hat{Y} &= \sum_{i=1}^N I_i \frac{1}{\pi_i} y_i  \ ,
\end{align}
where $I_i$ is one when a unit sampled and zero otherwise. Notice that
\begin{align} \label{sum3}
E\left( I_i \right) = \pi_i \ ,
\end{align}
and that the covariance of two units ($\mathrm{Cov}\left(I_i,I_j\right)$ for $i \neq j$) is not assumed to have any particular structure.

## Unbiasedness of Horvitz-Thompson estimator
Given:


\begin{enumerate}
\item $\frac{n}{N} = v$, for some $v \in (0,1)$,
\item $\frac{1}{\pi_i} > 0 \ \forall i \in \{1, ..., N\}$
\end{enumerate}

It is simple to see that the HT estimator of a sum is unbiased
\begin{align}\label{bias}
E(\hat{Y}) &= \sum_{i=1}^N  \frac{E\left( I_i \right)}{\pi_i}  y_i \\
           &= \sum_{i=1}^N y_i \\
           &= Y
\end{align}

## Decreasing variance of the Horvitz-Thompson estimator
Lemma: 
For random variables $A$ and $B$, without loss of generality, let the standard deviation of $A$ ($\sigma_A$) be less than the standard devation of $B$ ($\sigma_B$), so that $\sigma_A \geq \sigma_B$. Then, for some $k \in [0,1]$
\begin{align}\label{varb}
\sigma_B = k \sigma_A \label{assume1}
\end{align}

Claim:
\begin{align}
\mathrm{Var} \left(A + B\right)& \leq \sigma_A^2 \left(1 + 2 k + k^2\right) \ .
\end{align}

Proof:

By definition,
\begin{align}
\mathrm{Var} \left(A + B\right) &= \sigma_A^2 + \sigma_B^2 + 2 \sigma_{AB} \ ,  \label{mainvar2}
\end{align}
where $\sigma_{AB}$ is the covariance of $A$ and $B$. The formula that relates the correlation coefficient $\rho$ to the standard deviations and covaraince is
\begin{align}
\rho &= \frac{\sigma_{AB}}{\sigma_A \sigma_B} \\
\sigma_{AB} &= \rho \sigma_A \sigma_B
\end{align}
this is clearly maximized when $|\rho|=1$ and $\rho$ has the same sign as $\sigma_{AB}$. This means
\begin{align}
\sigma_{AB} & \leq \sigma_A \sigma_B
\end{align}
plugging that into \eqref{mainvar2} results in 
\begin{align}
\mathrm{Var} \left(A + B\right) & \leq  \sigma_A^2 + \sigma_B^2 + 2 \sigma_A \sigma_B \ ,
\end{align}
plugging in \eqref{assume1} shows 
\begin{align}
\mathrm{Var} \left(A + B\right) & \leq  \sigma_A^2 + k^2 \sigma_A^2 + 2 k\sigma_A^2 \ , \\
& \leq  \sigma_A^2  \left(1 + 2k +  k^2 \right) \ , \\
\end{align}

this completes the proof.

Now, consider the variance of the HT estimator

\begin{align}
\mathrm{Var}(\hat{Y}) &= \mathrm{Var} \left( \sum_{i=1}^N  I_i \frac{1}{\pi_i}  y_i \right)\ \label{vmat}
\end{align}

Consider building the variance unit by unit. First, for each unit
\begin{align}
\sigma_i^2 &= \frac{y_i^2}{\pi_i^2} \mathrm{Var}(I_i)
\end{align}
then, because $I$ has a Bernouli distribution its variance is always less than $\frac{1}{4}$ and
\begin{align}
\sigma_i^2 & \leq \frac{y_i^2}{4\pi_i^2}  \\
& \leq \mathrm{max} \frac{y_i^2}{4\pi_i^2} 
\end{align}
Note that the fact that $y_i \in \Re$ and $\pi_i > 0$ sufice to bound this max as a finite value.

For conveniance, then define
\begin{align}
q=\mathrm{max} \frac{y_i^2}{4\pi_i^2} \label{defq}
\end{align}

Then, consider building the variance estimate unit by unit.
\begin{align}
A_1&=\frac{I_1 y_1}{\pi_1} \\
B_1&=\frac{I_2 y_2}{\pi_2}
\end{align}
Then, by the lemma above,
\begin{align}
\mathrm{Var} \left( A + B \right) &\leq \sigma_A^2  \left(1 + 2k +  k^2 \right)
\end{align}
this is maximized when $k=1$, its largest possible value, so
\begin{align}
\mathrm{Var} \left( A + B \right) &\leq 4q 
\end{align}

Then, adding further units results in, for $m\geq2$
\begin{align}
A_m&=\sum_{i=1}^{m-1} \frac{I_i y_i}{\pi_i} \\
B_m&=\frac{I_m y_m}{\pi_m}
\end{align}
and, again, by the above lemma
\begin{align}
\mathrm{Var} \left( A + B \right) &\leq \sigma_{A_m}^2  \left(1 + 2k_m +  k_m^2 \right)\\
\end{align}
where $\sigma_{A_m}$ is a series defined by this eqation and knowing that $k_2 = 1$ and $\sigma_{A_2} = 4q$. this yeilds $k_3=1/4$ and $sigma_{A_3} = 4q (1+\frac{5}{16})$.







\begin{align}\label{vmat}
\mathrm{Var}(\hat{Y}) &= \sum_{i=1}^N \mathrm{Var} \left( I_i  \right) w_i^2  y_i^2 + \sum_{i=1}^{N-1} \sum_{j=i+1}^N w_i y_i  w_j y_j \mathrm{Cov}  \left( I_i, I_j  \right) \ .
\end{align}

Because the $I_i$ has a Bernouli distribution, $\mathrm{Var} \left( I_i  \right) \leq 1/4$.


To bound the covariance terms, consider building the variance of the sum elementwise. First, define
\begin{align}\label{def1}
A&=\sum_{i=1}^{m-1} y_i \\
B&=y_m
\end{align}
then 
\begin{align}
\mathrm{Var} \left[ \left(\sum_{i=1}^{m-1} y_i \right)  + y_m\right] &= \mathrm{Var} \left(A + B \right) \