classdef KMedianCoresetAlg < AbstractCoresetAlg
  %   Algorithm to compute the coreset for a given set of points in d
  %   dimentional space
  %   The algorithm works by taking a weighted set of points and finds a
  %   weighted coreset for these points.

  properties (Constant)
    linearInK = 1; % coresetType
    quadraticInK = 2; % coresetType
  end % (Constant)

  properties


    % Type = PointFunctionSet
    % Set of input points P
    P;

    % Type integer
    % k-median
    k;


    warningNegative;


    % this is the size of the coreset. It is used in sampling the
    % points with probabilities calculated after bicriteria algorithm.
    t;

    weightsFactor;

    % Coreset is a weighted set of points that is returned by the
    % algorithm
    % Type = PointFunctionSet
    coreset;

    report; % debugging and reporting information on the run

    coresetType;


    % added fields for keeping track of original indices
    sample_idx;
    num_bicriteria_points;

  end

  methods
    %constructor
    function obj=KMedianCoresetAlg()
      obj.bicriteriaAlg=AbstractBicriteriaAlg();
      obj.bicriteriaAlg.robustAlg.onSampleAlg = AllJsubspacesAlg;
    end
    % Compute the optimal cost of the set P
    optCost=computeOptCost(obj, P);

    % copmute the cost of the center Q to the set P
    cost=computeCost(obj, P, Q);

    % Compute the optimal solution of P
    opt=computeOpt(obj, P);


    % Merge two coresets C1, C2 to a new one.
    % C1 and C2 are PointfunctionSets. Compute the coreset on the merged
    % set of points. This is a merged Coreset.
    function C=mergedCoreset(obj, C1, C2)
      % Coresets are PointFunctionSets
      C1.merge(C2);
      % if isempty(obj.bicriteriaAlg)
      % C = obj.computeCoreset(C1);
      % else  C=obj.computeUsingBicriteria(SvdFunctionSet(C1,obj.j),obj.bicriteriaAlg);
      % end
      C = obj.computeCoreset(C1);
    end

    % Compute a coreset for input points. Coreset is a a weighted set of
    % points. This function should support weighted input point as it
    % will be applied recursively to find the coreset of coreset.
    % can pass isUpdate as second parameter
    function [C,sample_idx,num_bicriteria_points] = computeCoreset(obj, P)
      obj.P=P;
      %BClusterVec= obj.computeBicriteria();
      BClusterVec = ClusterVector;
      BClusterVec.initWithCenters( obj.P, PointFunctionSet( mean(obj.P.M.m)));
      if obj.coresetType == KMedianCoresetAlg.linearInK
        [Cmatrix, Cweights] = ...
          obj.computeCoresetLinearInK(BClusterVec);
      elseif obj.coresetType == KMedianCoresetAlg.quadraticInK
        [Cmatrix, Cweights]= ...
          obj.computeCoresetQuadraticInK(BClusterVec);
      else
        error ('wrong type of coreset');
      end % if
      obj.coreset = obj.matrixToFunctionSet(Cmatrix, Cweights);
      C=obj.coreset;

      % return sampled indices with the coreset
      sample_idx = obj.sample_idx;
      num_bicriteria_points = obj.num_bicriteria_points;

    end % computeCoreset

    function [Cmatrix, Cweights] = ...
        computeCoresetQuadraticInK(obj, BClusterVec)
      centers = BClusterVec.getCenterIndexes();
      clustersizes=BClusterVec.getClusterSize();
      additiveWeights = 8./clustersizes(centers);
      probs = 2*BClusterVec.getDistances./sum(BClusterVec.getDistances)+additiveWeights;
      [Sweights Sindices]=obj.samplePoints(probs);
      Cmatrix=obj.P.M.matrix(BClusterVec.pointIndexes(Sindices),:);
      Cweights=Sweights;
    end % computeCoresetQuadraticInK


    function [Cmatrix, Cweights] = ...
        computeCoresetLinearInK(obj, BClusterVec)

      % get number of points in each cluster.
      clustersizes=BClusterVec.getClusterSize();

      Bmatrix=BClusterVec.getCenters(1:BClusterVec.nClusters).M.m;

      probs = BClusterVec.getDistances;

      if sum(probs(:)) < 1e-10 % almost zero

        Smatrix = [];
        Sweights = [];
        SBw = 0;

        obj.sample_idx = [];
        obj.num_bicriteria_points = 0;

      else

        [Sweights Sindices]=obj.samplePoints(probs);

        obj.sample_idx = BClusterVec.pointIndexes(Sindices);
        obj.num_bicriteria_points = BClusterVec.nCenters;

        % Coreset is the union of the Sampled points and Bicriteria
        % points.
        Smatrix=obj.P.M.matrix(BClusterVec.pointIndexes(Sindices),:);

        % The inputs are weighted points. Multiply by the mp values
        % by weights before sampling.
        % Mp= (Mp).*(obj.P.W);
        % get the number of center from every item in Sindices
        Scenters=BClusterVec.getCenterIndexes(Sindices);



        % Get the points common in Sample S and the ith cluster. Find
        % the total weight of points for each such cluster.
        % SBw[i] contains the sum of weights of points common in S and
        % in the i-th cluster.
        b=BClusterVec.nCenters;

        if ~isempty(Scenters) && ~isempty(Sweights)
          SBw=accumarray(Scenters, Sweights', [b 1]);
        else
          SBw = 0;
        end

      end

      % Compute the weights for the bicriteria points.
      % Bicriteria center weights are equal to the number of points
      % in the cluster multiplied by (1+ epsilon) - ( the sum of
      % points that appera both in the cluster and the sample)
      Bweights=clustersizes-SBw;

      % TODO: remove points with zero weights

      % With high probability, Bweights should be positive, otherwise
      % we raise a warning flag
      neg=(Bweights<0);
      %if sum(neg)>0
      obj.warningNegative=true;
      Bweights(neg)=0;
      %else
      %    obj.warningNegative=false;
      %end

      % TODO Merge S and B. Nedd to do a union. Also can remove above
      % two statements.
      Cmatrix=[Smatrix; Bmatrix];
      Cweights=[Sweights; Bweights];
      %[Cmatrix iS iB]=union(Smatrix, Bmatrix, 'rows');

%       mydebug
%       disp(['size cmatrix = ' num2str(size(Cmatrix,1)) ',' num2str(size(Cmatrix,2)) ', ...
%         'size cweights = ' num2str(size(Cweights,1)) ',' num2str(size(Cweights,2))])

    end % computeCoresetLinearInK()



    % isUpdate is true or false
    function BClusterVec=computeBicriteria(obj)
      % Compute Bicriteria
      % Call the bicriteria algorithm. it returns you the set of
      % clusters which are used to compute the coreset.
      obj.report.bicriteriaTime=tic;


      % compute the bicriteria. This returns a Cluster vector
      % containing the centers and distances.
      BClusterVec=obj.bicriteriaAlg.compute(obj.P);
      obj.report.bicriteriaTime=toc(obj.report.bicriteriaTime);

      % Bicriteria running time.
      % Only for reporting
      obj.report.bicriteriaCost=sum(BClusterVec.getDistances);

      % Get the number of centers returned by the Bi-criteria
      % algorithm.
      obj.report.bicriteriaSize= BClusterVec.nCenters;
    end % computeBicriteria

    function [Sweights Sindices] = samplePoints(obj, probs)
      % Get set S by Non-Uniform Random Sampling the points. The set
      % S is used with the Bicriteria points to create the coreset.
      nrs=NonUniformSamplingAlg(obj.t);
      [Sweights Sindices]=nrs.sample(probs);
    end % samplePoints

  end

  methods (Static)
    % Todo: should be method of functionSet
    function C=matrixToFunctionSet(Cmatrix, Cweights)
      % Create a coreset Matrix
      CMatrix=Matrix(Cmatrix);
      CWeights=Matrix(Cweights);

      % Create a pointFunction Set for the Coreset
      C=PointFunctionSet(CMatrix);
      C.W=CWeights;
    end % matrixToFunctionSet
  end
end

function [sqdists]=sqDistances(q, points)

[n,p] = size(points);
dists = zeros(n,2);
sqdists = zeros(n,1);

dists(:,1) = (points(:,1) - q(1,1)).^2;
dists(:,2) = (points(:,2) - q(1,2)).^2;
sqdists=sum(dists, 2);
end % function [sumSqDists]

