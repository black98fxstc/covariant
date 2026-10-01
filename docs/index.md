---
layout: home
title: Leonard
---

Leonard is a desktop tool for the analysis of flow cytometry datasets using
methods from computational and differential geometry.

## Covariant Statistics

Covariant Statistics is the geometric foundation behind Leonard’s clustering
and projection-pursuit workflow. It treats a sample as a high-dimensional
distribution and analyzes how that distribution changes under scale-preserving
transformations.

At a high level, Leonard uses these ideas to:

- factor an *n*-dimensional sample distribution into many one-dimensional views,
- identify regions where the Laplacian is positive and the data behave like a
  smooth cluster or phenotype,
- separate ambiguous boundary regions where the tendency is toward dispersion,
- and compute invariant quantities that remain meaningful across different data
  scales and monotonic transformations.

For a fuller explanation, see [Covariant Analysis]({{ '/covariant/' | relative_url }}).

## Analysis Methods

- **[Exhaustive Projection Pursuit]({{ '/methods/' | relative_url }}#exhaustive-projection-pursuit)**: Examines all two-dimensional projections to find the best split into two new populations, the simplest and most interpretable case.
- **[Laplacian Clustering]({{ '/methods/' | relative_url }}#laplacian-clustering)**: An *n*-dimensional clustering method based on ideas from differential geometry (supporting 2, 3, and 4 dimensional cases).
- **[Covariant Statistics]({{ '/methods/' | relative_url }}#covariant-statistics)**: Statistics invariant under monotonic transformations (independent of data scale), offering a quantitative morphology layer on top of clustering.

## Workflow

Leonard can be run from the command line on a variety of data file formats but primarily operates on FlowJo workspaces (`.wsp`), enabling the analysis of FlowJo gated subpopulations. The investigation typically follows a simple flow: open a workspace, choose the sample and population of interest, select the analytical dimensions, and generate reports.

## Explore Leonard

- [About Leonard]({{ '/about/' | relative_url }}): Background, history, and project dedication.
- [Analysis Methods]({{ '/methods/' | relative_url }}): Detailed explanations of computational and differential geometry methods.
- [Documentation]({{ '/documentation/' | relative_url }}): System requirements, workflows, and usage guide.
- [Example Reports]({{ '/examples/' | relative_url }}): Overview and links to example analysis reports.

## Requirements

macOS and Windows are supported. Computers suitable for running FlowJo should be fine.

## Downloads & Installation

- **Downloads**: Get the latest alpha builds from [GitHub Releases](https://github.com/black98fxstc/covariant/releases/latest).
- **Installation**: See the [installation instructions](https://github.com/black98fxstc/covariant#readme).
- **Source Code**: [black98fxstc/covariant repository](https://github.com/black98fxstc/covariant).
- **Security**: [Security policy](https://github.com/black98fxstc/covariant/security/policy).

## References

- Moore, W. A., Meehan, S. W., Meehan, C., Parks, D. R., Walther, G., & Herzenberg, L. A. (2025).  
  *Automatic phenotyping using exhaustive projection pursuit*.  
  [Communications Biology](https://doi.org/10.1038/s42003-025-08581-z).  
  PMID: [40797028](https://pubmed.ncbi.nlm.nih.gov/40797028/) · PMCID: PMC12343891

- Wayne A. Moore.  
  *Covariant Statistics*.  
  [Current Draft](https://drive.google.com/file/d/1Sq5W18-JzMJAaWGBl-vYJojl5mElBpw4/view?usp=drive_link)
