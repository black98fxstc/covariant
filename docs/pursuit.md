---
layout: page
title: Exhaustive Projection Pursuit
permalink: /pursuit/
---

### Exclude uniformly stained dimensions

First 1-dimensional distributions of the sample are compared to a uniformly stained distribution using the [Kullback-Leibler Divergence](https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence) and only those showing significant structure are considered.

### Sample density estimator

For each pair of dimensions the events are projected down onto the corresponding plane and the [weights]{https://doi.org/10.1155/2009/686759}, a more sophisticated form of histogram, are calculated. The weights are then smoothed with a [Gaussian filter](https://en.wikipedia.org/wiki/Gaussian_filter) to form a density estimator.

### Find all possible candidate separations

The density estimator is subject to modal clustering, which produces a plane graph of vertices edges and faces, one face for each mode (peak) in the sample distribution and the edges are the boundaries between them. We always want to divide into only two parts, because other dimensions than the two under consideration may provide a superior split to the additional edge in this pair. If there are more than two clusters, every possible simple subgraph must be tried.

### [Density based merging](https://doi.org/10.1155/2009/686759)

Due to limited sample size, modal clusters may be found that are the results of random sampling errors. To guard against this, a check is made to ensure that the dip between two peaks peaks is statistically significant. If it is not, the edge is removed from the cluster graph.

### The fewer events near the boundary the better

The splits are scored according to the number of events close to the boundary. When all the possible splits of all the pairs are complete, the best candidate found (if any) is chosen to divide the sample and the process starts over on each part.