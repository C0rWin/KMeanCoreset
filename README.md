k-Means for Streaming and Distributed Big Sparse Data
===

Introduction
---

In this repository we provide Matlab implementation of coreset algorithms, used
for evaluation in following [paper](https://arxiv.org/pdf/1511.08990.pdf):

> **k-Means for Streaming and Distributed Big Sparse Data.**
> *Artem Barger, and Dan Feldman.*
> Proceedings of the 2016 SIAM International Conference on Data
> Mining. Society for Industrial and Applied Mathematics, 2016.

Algorithms
---

We provide implemenation for three algorithms used in the paper above:

1. **Uniform coreset**
2. **Non-uniform corest (sensitivity based)**
3. **Our algorithm (Deterministic coreset construction)**

API
---

Coreset algorithms provides two very basic API's:

* Given set of points P from `R^d`:

`computeCoreset(P`) - compress points P into coreset the weighted set C

* Given two coresets C1 and C2:

`mergedCoreset(C1, C2)` - merges two coresets into a new one C'.

Detailed Usage
---

1. **Matrix** - `Matrix.m`

Matrix abstraction which encapsulates a set of `P` points of `n` in `R^d` in
matrix of size `n-by-d`.

2. **PointFunctionSet** - `PointFunctionSet.m`

Class which represent weighted point set, in terms of function which maps point
into real value (weight).

3. **Uniform coreset** - `uniformCoreset.m`

Implementation of uniform "naive" coreset sampling, with following API's

```matlab
coresetSize = 100;

algorithm = uniformCorest(coresetSize);

% Compute coreset of n points from R^d
coreset1 = algorithm.computeCoreset(P1);
coreset2 = algorithm.computeCoreset(P1);

% Merge two coresets into new one
coreset = algorithm.mergedCoreset(coreset1, coreset2);
```

4. **Non uniform coreset (sensitivity based)** - `KMedianCoresetAlg.m`

Implementation of sensitivity based coreset sampling - the non uniform, for
both `k-means` and `k-median` algorithms. Algorithms uses bicriteria
approximation for sensitivity computation.

```matlab
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

5. **Determenistic coreset algorithm** - `kmeansCoreset.m`

Implementation of our determenistic kmeans coreset algorithm of size k^O(1).

```matlab
coresetSize = 100;
% Number of kmeans++ iterations to execute for coreset construction
maxIter = 10;

algorithm = kmeansCorest(coresetSize, maxIter);

% Compute coreset of n points from R^d
coreset1 = algorithm.computeCoreset(P1);
coreset2 = algorithm.computeCoreset(P1);

% Merge two coresets into new one
coreset = algorithm.mergedCoreset(coreset1, coreset2);

```

6. *Coresets streaming* - `Stream.m`

Stream encapsulates logic of building coreset merge-and-reduce tree and allow
streaming computation of the coreset for continiously arriving points.

Example of streaming:

```matlab

% We will use kmeans coreset as an example here:

coresetSize = 100;
maxIter = 10;

algorithm = kmeansCorest(coresetSize, maxIter);

stream = Stream();
stream.coresetAlg = algorithm;
% Define coreset merge-and-reduce tree leaf size
stream.leafSize = 100;

% For each new batch of streaming points do:
stream.addPointSet(P)

% In order to get final coreset results, which is union of entire tree into
final single coreset of given size.

result = stream.getUnifiedCorest()

```

7. **Evaluation**

Evaluation of results done with following code:

```matlab
% Given kmeans results computed with coreset, we compute the value of energy
% function, which is sum of squared distances to the kmeans centers and compare
% with kmeans++ approximated solution for original dataset

coresetSize = 100;
% Number of kmeans++ iterations to execute for coreset construction
maxIter = 10;

algorithm = kmeansCorest(coresetSize, maxIter);

% Compute coreset of n points from R^d
coreset = algorithm.computeCoreset(P);

% coreset results is of type PointFunctionSet, it has matrix which represents
% coreset points and well matrix which actually holds weights of the points
% computed by algorithm.

[part, centers, ~, ~] = Ckmeans(coreset.M.m, k, coreset.W.m, 'distance', 'sqeuclidean', ...
        'maxiter', 100, 'emptyaction', 'singleton', 'display', 'off', 'onlinephase', 'off');

dist = zeros(size(P, 1), k);
for c=1:k
    dist(:,c) = sum(bsxfun(@minus, P, centers(c,:)).^2, 2);
end
energy = sum(min(dist, [], 2));

error = energy / optEnergy - 1
```

Feedback
---

If you'd like to use this implementation, please reference original paper, any
feedback send to Artem Barger (artem@bargr.net)

License
---

The software is released under the MIT License as detailed in kmeans.pyx.

