---
layout: page
title: Covariant Analysis
permalink: /covariant/
---
### Factoring Probability

The first step is to construct a sample density estimator as described for [Exhaustive Projection Pursuit]({{ '/pursuit/' | relative_url }}#sample-density-estimator). The n-dimensional sample distribution function is factored into a large number of 1-dimensional distributions, which are statistically independent. 

### Natural Parameters

These 1-dimensional functions must obey a differential equation, which is used to find the first and second derivatives of the logarithm of these density functions. The results are known in statistics as the natural parameters and in geometry as the Christoffel symbols.

### Laplacian Clustering

The Laplacian is the sum of some of these second derivatives. When the Laplacian of the distribution function of a population is positive, it has a smooth distribution with a strong central tendency, which is to say, it looks like what we think of as a cluster or phenotype. The events in a contiguous region with a consistent positive Laplacian define a Laplacian cluster. The clusters are surrounded and separated by regions where the Lagrangian is negative and the tendency is to dispersion not concentration. Events in this region are ambiguous and cannot be assigned to any specific cluster.

### Covariant Statistics

There exist invariant quantities that do not depend on the scales used. For a common example, you would get the same answer if you used either linear or logarithmic scales. The first example was found by Gauss in the nineteenth century, a result he called Theorema Eggregium (Remarkable Result) it was so surprising.

The Ricci and Einstein tensors (yes that Einstein) are also sums of different subsets of the second derivatives. In this case the invariant values are sums of these over the events in each Laplacian population. For any normal distribution, all these vanish, so in a way they measure deviation from normality. In particular, they measure the precision with which you can predict one variable, based on knowing the value of another. When the invariants are not zero, the sample does not have a single fixed mean but rather means distributed along some curve.

It is substantially faster to compute the Laplacian alone, therefore it is implemented as a separate step for convenience when only the clusters are wanted.
