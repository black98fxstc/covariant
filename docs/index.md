# Leonard

Leonard is a desktop tool for the analysis of flow cytometry datasets using
methods from computational and differential geometry, inluding

- Exhaustive Projection Pursuit
- Laplacian Clustering
- Covariant Statistics

It can be run from the command line on a variety of data file formats but
primarily will be run on FlowJo workspaces (.wsp) and then allows the analysis
of FlowJo gated subpopulations. The investigator selects a sample, a set of 
dimensions and then one or more previously identified populations.

Exhaustive Projection Pursuit examines all two dimensional projections to find
the best split into two new populations, then starts the process on each subpopulation 
until no suitable split can be found. It is suitable for any number of dimensions,
and it is possible to sort on the resulting gating tree.

Laplacian Clustering is an n-dimensional method based on ideas from differential
geometry. 2, 3, and 4 dimensional analysis is supported and while the subpopulations
may be further analyzed in software, existing hardware does not support sorting.

Covariant Statistics are statistics that are invariant under monotonic transformations,
i.e., independant of the data scale used. They offer a quantitative morphology of
Laplacian cluster shapes that measure how strongly the values of different dimensions 
are linked.

## Requirements
Mac and Windows are supported. Computers suitable for running FlowJo should be fine.

## Downloads
Get the latest alpha builds from [GitHub Releases](https://github.com/black98fxstc/covariant/releases/latest)

## Installation
[Installation instructions](https://github.com/black98fxstc/covariant/blob/main/README.md)

## Source code
[Source code repository](https://github.com/black98fxstc/covariant)

## Security
[Security policy](https://github.com/black98fxstc/covariant/security/policy)

## References

- Moore, W. A., Meehan, S. W., Meehan, C., Parks, D. R., Walther, G., & Herzenberg, L. A. (2025).  
  *Automatic phenotyping using exhaustive projection pursuit*.  
  [Communications Biology](https://doi.org/10.1038/s42003-025-08581-z).  
  PMID: [40797028](https://pubmed.ncbi.nlm.nih.gov/40797028/) · PMCID: PMC12343891

- Wayne A. Moore
  *Covariant Statistics*. 
  [Current Draft](https://drive.google.com/file/d/1Sq5W18-JzMJAaWGBl-vYJojl5mElBpw4/view?usp=drive_link)

## About
[About Leonard](https://github.com/black98fxstc/covariant/ABOUT.md)
