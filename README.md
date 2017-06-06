k-Means for Streaming and Distributed Big Sparse Data
===

Introduction
---

In this repository we provide Matlab implementation of coreset algorithms, used
for evaluation in following paper:

> **k-Means for Streaming and Distributed Big Sparse Data.**
> *Artem Barger, and Dan Feldman.*
> Proceedings of the 2016 SIAM International Conference on Data
> Mining. Society for Industrial and Applied Mathematics, 2016.

Algorithms
---

We provide implemenation for three algorithms used in the paper above:

1. Uniform coreset
2. Non-uniform corest (sensitivity based)
3. Our algorithm (Deterministic coreset construction)

Detailed Usage / API
---

1. Matrix - `Matrix.m`

Matrix abstraction which encapsulates a set of `P` points of `n` in `R^d` in
matrix of size `n-by-d`.

2. PointFunctionSet - `PointFunctionSet.m`

Class which represent weighted point set, in terms of function which maps point
into real value (weight).

3. Uniform coreset - `uniformCoreset.m`

Implementation of uniform "naive" coreset sampling, with following API's

```
coresetSize = 100;

algorithm = uniformCorest(coresetSize);

% Compute coreset of n points from R^d
coreset1 = algorithm.computeCoreset(P1);
coreset2 = algorithm.computeCoreset(P1);

% Merge two coresets into new one
coreset = algorithm.mergedCoreset(coreset1, coreset2);
```

4. Non uniform coreset (sensitivity based) - `KMedianCoresetAlg.m`

Implementation of sensitivity based coreset sampling - the non uniform, for
both `k-means` and `k-median` algorithms. Algorithms uses bicriteria
approximation for sensitivity computation.

```
coresetSize = 100

algorithm = KMedianCoresetAlg();

% For k-median use KMedianCoresetAlg.linearInK
algorithm.coresetType = KMedianCoresetAlg.quadraticInK;

% Setup bicriteria algorithm parameters, basically configure
% alpha and beta parameters of the approximation.

algorithm.bicriteriaAlg.robustAlg.beta = beta;
algorithm.bicriteriaAlg.robustAlg.partitionFraction = alpha;
algorithm.bicriteriaAlg.robustAlg.costMethod = ClusterVector.sumSqDistanceCost;

% Compute coreset of n points from R^d
coreset1 = algorithm.computeCoreset(P1);
coreset2 = algorithm.computeCoreset(P1);

% Merge two coresets into new one
coreset = algorithm.mergedCoreset(coreset1, coreset2);
```

Feedback
---

If you'd like to use this implementation, please reference original paper, any
feedback send to Artem Barger (artem@bargr.net)

License
---

The software is released under the MIT License as detailed in kmeans.pyx.

