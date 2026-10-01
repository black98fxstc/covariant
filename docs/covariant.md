---
layout: page
title: Covariant Analysis
permalink: /covariant/
---

The first step is to factor the n-dimensional sample distribution function into 
a large number of 1-dimensional distribution functions. The first and second 
derivatives of these functions are known in statistics as the natural parameters 
and in geometry as the Christoffel symbols.

The Laplacian is the sum of some of these second derivatives. When the Laplacian 
of a distribution function of a population is positive, it is a smooth distribution 
with a strong centeral tendency, which is to say, that it looks like what we think 
of as a cluster or phenotype. The events in a contiguous region with a consisten 
positive Laplacian define a Laplacian cluster. The clusters are surrounded and 
separated by areas where the Lagrangian is negative and the tendency is to dispersion 
not concentration. Events in this reageon are ambiguous and cannot be assighed to 
any specific cluster.

There exist invariant quantities that do not depend on the scale chosen, i.e., you 
get the same answer if you plot the data linearly that you do if you plot it 
logarithmically. The first example was found by Gauss in the nineteeth century a 
result he called Theorema Eggregium (Remarkable Result) it was so surprising.

The Ricci and Einstein tensors (yes that Einstein) are also sums of different subsets 
of the second derivatives. In this case the invariant values are sums of these over 
the events in each Laplacian population. For any normal distribution all these vanish, 
so in a way they measure deviation from normality. In particular, they measure the 
precision with which you can predict one variable, based on knowing the value of another.

It is substantially faster to compute the Laplacian alone, therefore it is implemented 
as a separate step for convenience when only the clusters are wanted.