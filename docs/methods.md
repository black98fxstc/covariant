---
layout: page
title: Methods
permalink: /methods/
---

Leonard brings methods from computational and differential geometry to the analysis of flow cytometry data.

## Exhaustive Projection Pursuit

Exhaustive Projection Pursuit examines all two-dimensional projections of a multi-dimensional dataset to identify optimal splits into distinct subpopulations. The process recurses on each subpopulation until no suitable split can be found.

- **Dimensionality**: Suitable for any number of dimensions.
- **Physical Sorting**: Gating trees produced by projection pursuit can be used directly for sorting on compatible cytometer hardware.
- **Reference**: Moore, W. A. et al. (2025). *Automatic phenotyping using exhaustive projection pursuit*. [Communications Biology](https://doi.org/10.1038/s42003-025-08581-z).

## Laplacian Clustering

Laplacian Clustering is an *n*-dimensional clustering method rooted in differential geometry.

- **Dimensionality**: Supports 2, 3, and 4-dimensional analysis.
- **Analysis**: Subpopulations may be further explored and analyzed in software.

## Covariant Statistics

Covariant Statistics are statistical measures invariant under monotonic transformations (independent of the data scale used).

- **Morphology**: Offers a quantitative morphology of Laplacian cluster shapes.
- **Linkage**: Measures how strongly values of different dimensions are linked.
- **Reference**: Wayne A. Moore. *Covariant Statistics*. [Current Draft](https://drive.google.com/file/d/1Sq5W18-JzMJAaWGBl-vYJojl5mElBpw4/view?usp=drive_link).

## Next Steps

- Review the [Documentation]({{ '/documentation/' | relative_url }}) for workflow details.
- View [Example Reports]({{ '/examples/' | relative_url }}).
