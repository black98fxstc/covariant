---
layout: page
title: Exhaustive Projection Pursuit
permalink: /covariant/
---

First 1-dimensional distributions of the sample are compared to a uniformly stained distribution 
using the Kuhlbach-Leibler Divergence and only those showing some structure are considered.

For each pair of dimensions the events are projected down onto the corresponding plane and 
the weights, a more sophisticated form of histogram, are calculated. The weights are then 
smoothed with a gaussian kernel to form a density estimator.

The density estimator is subject to modal clustering, which produces a plane graph of verticies,
edges and faces, one face for each mode (peak) in the sample distribution. We always want to 
divide into only two parts, because other dimensions than the two under consideration may provide 
a superior split to the additional edge in this pair. If there are more than two clusters, every 
possible simple subgraph must be tried.

The splits are scored accourding to the number of events close to the boundary. When all the possible 
splits of all the pairs are complete the best candidate found (if any) is chosen to divide the 
sample and the process starts over on each part.