classdef ClusterVector < HandleObject
%--------------------------------------------------------------------------
%This class is a general implementation of the set of clusters. Although it
%is designed for the bicriteria clustering algorithm at the very beginning,
%in fact it is a general representation of any set of clusters by any
%clustering algorithms. In order to make the class self-contained, the
%major attributes we need to maintain for any set of clusters are:
%   - The set of points to be clustered
%       - Corresponds to the member variable "F";
%
%   - The set of centers in the clusters (we use the indexes of the centers to name the clusters)
%       - Corresponds to the member variable "centers";
%       - Member variable "centerCounter" records the number of centers in
%       the clusters (or the number of non-zero elements in centers);
%
%   - The indexes of points which are already clustered (in some applications we may not want to cluster the whole set)
%       - Corresponds to the member variable "pointIndexes";
%       - Member variable "pointCounter" records the number of clustered
%       points (or the number of non-zero elements in pointIndexes);
%
%   - The index of the cluster(center) a point belongs to
%       - Corresponds to the member variable "centerIndexes";
%
%   - The distance of a point to its cluster center 
%       - Corresponds to the member variable "distances";
%
%Here is an example:
%  
%Initially, suppose F contains 10 points. Naturally, the indexes of the 10
%points are from 1 to 10. In the following I shall use the index of every
%point as its name. After running some clustering algorithm (not restricted
%to the bicriteria algorithm), suppose we have the following centers:
%   centers: point 2, point 5, point 7, point 10. 
%
% For each center, it defines a cluster with some points in F associated
% with the center. Suppose we have the following clusters ("dist" means the
% distance from the point to its center):
%   Cluster 1: center: point 2;
%              points in the cluster: point 2(dist: 0), point 3(dist: 0.4),
%              point 6(dist: 0.9), point 8(dist: 1.0);
% 
% 
%   Cluster 2: center: point 5; 
%              points in the cluster: point 1(dist: 0.3), point 4(dist:
%              0.7), point 5(dist: 0);
% 
%   
%   Cluster 3: center: point 7;
%              points in the cluster: point 7(dist: 0), point 9(dist: 0.4);
% 
%   
%   Cluster 4: center: point 10;
%              points in the cluster: point 10(dist: 0);
% 
% 
%For this case, the member variables in the class looks as follows: 
%   pointIndexes: [2, 3, 6, 8, 1, 4, 5, 7, 9, 10]
%       Short explanation: we rearrange the points in F according to
%       clusters and put them into pointIndexes from Cluster 1 to Cluster
%       4. Every value in pointIndexes corresponds to the original index of
%       the point in F. (Recall that we use the indexes in F as the name
%       for every point)
% 
%   distances: [0, 0.4, 0.9, 1.0, 0.3, 0.7, 0, 0, 0.4, 0]
%       Short explanation: the distances vector holds a one-to-one
%       correspondence with the pointIndexes vector, i.e., distances[i] (1
%       <= i <= 10) is the distance between point pointIndexes[i] to its
%       corresponding center.
% 
%   centerIndexes: [1, 1, 1, 1, 2, 2, 2, 3, 3, 4]
%       Short explanation: the centerIndexes vector also holds a one-to-one
%       correspondence with the pointIndexes vector, i.e., centerIndexes[i]
%       (1 <= i <= 10) is the index of center/cluster the point
%       pointIndexes[i] belongs to. Note that here we change the indexes
%       for centers. For example, although the center of Cluster 1 is point
%       2, in the record of centerIndexes, I take "point 2" as "center 1".
%       In order to be able to fetch the "real index" of a center, we
%       provide the following "centers" member variable.
% 
%   centers: centers[i] (1 \leq i \leq 4) refers to the coordinates of
%       "center i". For example, center[2] refers to the coordinates of
%       "center 2", which is in fact the coordinates of point 5.
% 
%   F: F is still the intact set of points. F[i] (1 \leq i \leq 10) refers
%       to the coordinates of point i.
% % 
%   centerCounter: since we have 4 centers corresponding to 4 clusters,
%       centerCounter = 4.
%
%--------------------------------------------------------------------------
    properties (Constant)
        % cost function for approximation
        maxDistanceCost = 1;
        sumDistanceCost = 2;
        sumSqDistanceCost = 3;
    end

    properties
        distances; % distances[i]: distance between pointIndexes[i] to its corresponding center.
    end

    properties (SetAccess = protected)
        % To understand better of the meaning of all the following member
        % variables, please refer to the above example.        
        
        pointIndexes; % Record the indexes of points in the original point set F which have been clustered.

        centerIndexes; % centerIndexes[i]: the center index/cluster index of point pointIndexes[i]. 
                       % Note that the center index is different from the point index in F. 
                       % We shall use the following member variable
                       % "centers" to connect the two.

        % Type =  MATLAB matrix
        centers; 
    
        % Type = PointFunctionSet
        F = PointFunctionSet; % F.M.m[i]: the coordinates of "point i";
           % F.W.m[i]: the weight of "point i". 
    end % private properties
    
    properties (Dependent)
         % Type = integer
         % Number of points to be clustered
         nClusters;
         nPoints;
         
         % Type = integer
         % Point dimension
         d;   
    end
     
    methods
        function set.distances(obj, dists)
            if size(dists,2)~=1
                % put a breakpoint here
                error('a');
            end
            obj.distances=dists;
        end

        function clustersizes=getClusterSize(obj)
            b=obj.nCenters;
            clustersizes=histc(obj.getCenterIndexes,1:b);
        end
        function weighted=getWeightedCenters(obj)
            weighted=obj.centers;
            weighted.W.m=obj.getClusterSize()';
            weighted.W.m=obj.normalizeWeights(weighted.W.m);
        end
        function weight=normalizeWeights(obj, weights)
            for i=1:length(weights)
                weights(i,1)=log(weights(i,1)+10);
            end
            weight = weights;
            %sumWeights = sum(weights)
            %avgWeights = sumWeights/length(weights);
            %weight = weights/avgWeights;
        end
        function s=asString(obj)
            s='testMe';
        end
        function result=isInfDist(obj)
            if sum(obj.distances==inf)>0
                result=true;
            else
                result=false;
            end
        end
        function cluster = epsilonGrid(obj, eps)
            s=[];
            isFirst=true;
            for i=1:obj.centers.size
                %[~, idx] = obj.getClusterIndexes(i);
                [temp, idx] = obj.getClusterIndexes(i);
                if ~isempty(idx)
                    [center, clusterIndexes] = ...
                        obj.getClusterIndexes(i);
                    G=obj.subsetPoints(clusterIndexes);
                    subs=center.epsilonGrid(eps, G);
                    if isFirst
                        s=subs;
                        isFirst=false;
                    else
                        s.merge(subs);                    
                    end; % if
                end % if
            end % for
            cluster = ClusterVector;
            cluster.initWithCenters(obj.F, s);
        end
        
        function clusterVec = subsetApply(obj, f)
            subsetIdx=[];
            for i=1:obj.centers.size
                [center, idx] = obj.getClusterIndexes(i);
                if ~isempty(idx)
                    G=obj.subsetPoints(idx);
                    chosenIdx = f(center, G, i);
                    chosenIdx =  idx(chosenIdx);
                    % should use pointIndexes instead
                    subsetIdx=[ subsetIdx; chosenIdx];
                end
            end
            clusterVec=obj.subset(subsetIdx);
        end
        
        % set point functions on subset of (local) indexes
        function setPointsSubset(obj, localIndexes, G)
            if G.M.d~=obj.F.M.d || G.size~=length(localIndexes)
                error('size do not match');
            end
            indexes=obj.pointIndexes(localIndexes);
            obj.F.setSubset(indexes, G);
        end
        % new copy of function set
        function P=subsetPoints(obj, indexes)
            idx=obj.pointIndexes(indexes);
            P=obj.F.subset(idx);
        end
        
        function clusterVec = apply(obj, f)
            clusterVec=ClusterVector();
            clusterVec.copy(obj);
            isFirst=true;
            for i=1:obj.centers.size
                [center idx] = obj.getClusterIndexes(i);
                if ~isempty(idx)
                    newPoints = f(center, obj.subsetPoints(idx));
                    if isFirst
                        clusterVec.F.M.m=NaN(clusterVec.F.n, newPoints.M.d);
                        clusterVec.F.W.m=NaN(clusterVec.F.n, 1);
                        isFirst=false;
                    end
                    clusterVec.setPointsSubset(idx, newPoints);
                end
            end % for
        end
        
        % Replace F with only relevant functions in the cluster
        %function compress(obj)
            %G=obj.F.copy();
            %G.M = obj.F.M.getRows(obj.pointIndexes);
            %G.W = obj.F.W.getRows(obj.pointIndexes);
            %initWithFunctionSet(obj, G, 1:length(obj.pointIndexes), obj.distances, obj.centers, obj.centerIndexes)
        %end
        function copy(obj, other)
           obj.initWithFunctionSet(other.F, other.pointIndexes, other.distances, other.centers, other.centerIndexes)
        end
        function result=nCenters(obj)
            result=obj.centers.size;
        end
        function result=getFunctionSet(obj, isSorted)
            if ~exist('isSorted', 'var')
                isSorted = false;
            end
            if isSorted
                idxes=sort(obj.pointIndexes);
            else
                idxes=obj.pointIndexes;
            end
            result = obj.F.subset(idxes);
        end %get the entire Clustered functionSet
        
        function result=getSumWeights(obj)
            G=obj.F.subset(obj.pointIndexes); 
            result=G.sumW;
        end %get sum of weights of a set of points in F w.r.t pointIndexes
        
        function result=getDistances(obj)
            result=obj.distances;
        end % get the distances vector 
        function result=getWeightedDistances(obj)
            weights = obj.F.W.m(obj.pointIndexes); 
            result=obj.distances.*weights;
        end % get the distances vector 

        function result=getPointCluster(obj,i)
            iclusterIndexes = (obj.centerIndexes == i);
            clusterIndexes = obj.pointIndexes(iclusterIndexes);
            result=obj.F.subset(clusterIndexes);
        end

        function result=getCenters(obj,indexes)
            if nargin >1
                result = obj.centers.subset(indexes);
            else
                result=obj.centers.subset(1:obj.centers.size);
            end
                
        end % get the centers
%         
%         % Initialize a new cluster vector
%         % It is only called in the cluster vector where you want to store
%         % the final results.
%         % For temporary cluster vectors, you don't need to call it as it
%         % costs much to clear and initialize
%         function init(obj, centerSize, intNandD)
%             if (centerSize > 0)
% %                 obj.centers = zeros(centerSize, obj.d);
%             end
%             
%             if (intNandD == true)
%                 obj.clear();
%                 % Since we can know in prior the size of all the following
%                 % vectors. It is faster to allocate enough space initially.
%                 obj.distances = zeros(obj.n, 1);      
%                 obj.centerIndexes = zeros(obj.n, 1);
%                 obj.pointIndexes = zeros(obj.n,1);
%             end
%         end
% 

        % init cluster with points and their centers
        function initWithCenters(obj, F, centers)
            obj.F = F;
            obj.centers = centers ; 
            obj.pointIndexes = 1:F.M.n;
            obj.updateNearestCenters();
        end

        % initial the cluster vector with a given set of clusters with
        % centers
        function initWithFunctionSet(obj, F, pointIndexes, distances, centers, centerIndexes)
            [r,c] = size( distances);            
            if ( c ~= 1) 
                r
                c
                error( 'Distances has to be be column vector');
            end
            obj.F=F.copy();
            obj.pointIndexes = pointIndexes;
            obj.distances=distances;
            obj.centerIndexes=centerIndexes;
            obj.centers=centers;
            
            % ensure centerindexes is a column vector
            if size(obj.centerIndexes,1) == 1
              %mydebug
              obj.centerIndexes = obj.centerIndexes';
            end
            
        end % initWithFunctionSet
        
        % Return clusterVec of the smallest partitionFraction of the
        % current
        function result = smallest(obj, partitionFraction)           
            % Calculate the index m such that the sum of weights of index
            % mIndexes is at least the "partitionFraction" of the total weights.
            dists=obj.distances;
            G=obj.F.subset(obj.pointIndexes); 
            [mIndexes os]= ...
                OrderStatisticBySortAlg.compute(dists, G.sumW * partitionFraction, G.W.matrix);
            if os==inf
                error('inf os');
            end
            result=obj.subset(mIndexes);
        end

%         % Append the PointFunctionSet to cluster.
%         % No duplicate check.
%         function obj=unionPoints(obj, P)
%             c=ClusterVector();
%             c.initWithCenters(P, obj.centers);
% 
%             n=obj.F.n;
%             obj.F.merge(P);
%             newindexes=n+1:obj.F.n;
%             obj.union(c, newindexes);
%         end

        % merge disjoint clusterVec
        function merge(obj, otherVec)
            otherN = otherVec.F.size();
            currentIndexes = (1:otherN)+obj.F.size;
            obj.F.merge(otherVec.F);
            obj.union(otherVec, currentIndexes);
        end
        
        % Append the specific cluster vector to the current cluster vector.
        % clusterVec.F should be a subset of obj.F:
        % clusterVec.F=obj.F.subset(currentIndexes);
        % but may use different centers.
        % No check for duplicate points.
        % if sameCenters==true, don't use new centers
        function obj=union(obj, clusterVec, currentIndexes, sameCenters)
            if isempty(obj.centers) % current cluster is empty
                obj.initWithFunctionSet(clusterVec.F, ...
                    clusterVec.pointIndexes, clusterVec.distances, ...
                    clusterVec.centers, clusterVec.centerIndexes)
            else
                if ~exist('sameCenters', 'var')
                    sameCenters=false;
                end
                if ~sameCenters
                    ocC = obj.centers.size;
                    % Append "centers" to the "centers" of the object
    %                clusterVec.centers.M.matrix
                    obj.centers.merge(clusterVec.centers);
                else
                    ocC=0;
                end
                
                % Use intermediate variables to speedup
                opC = length(obj.pointIndexes);
                cvpC = length(clusterVec.pointIndexes);
                
                % Append "distances" to the "distances" vector of the object
                obj.distances(opC + 1: opC + cvpC,1) = clusterVec.distances;
                
                % Append "centerIndexes" to the "centerIndexes" of the object
                obj.centerIndexes(opC + 1: opC + cvpC,1) = clusterVec.centerIndexes + ocC;
                
                % Append "pointIndexes" to the "pointIndexes" of the object
                try
                    obj.pointIndexes(1,opC + 1: opC + cvpC) = currentIndexes(clusterVec.pointIndexes);
                catch
                    error('df');
                end
            end
        end % method union
        
        % Recompute the point distances to the nearest centers
        function dists = updateNearestCenters(obj)
            G = obj.F.subset(obj.pointIndexes);
            [obj.distances, obj.centerIndexes] = G.Eval(obj.centers);
            dists = Matrix(obj.distances);
        end
        
        function setDistances(obj, d)
            obj.distances = d;
        end
        
        % Return a specific cluster i
        % center: the center of cluster i
        % clusterIndexes: indexes of points in the original F that belongs
        % to cluster i
        function [center localIndexes] = getClusterIndexes(obj, i)
              if (i > obj.centers.size)
                  disp('Cluster number exceeded!');
              else
                 iclusterIndexes = (obj.centerIndexes == i);
                 localIndexes=find(iclusterIndexes);
                 center = obj.centers.getAt(i);
                 %clusterIndexes = obj.pointIndexes(iclusterIndexes);
              end
        end
          
        %Return a specific cluster (subset of this obj)
        %currently works only for centers of type cell
        function [subClusterVec subCenter] = getCluster(obj,index)
            % use getAt only if centers is not array  
               subCenter = obj.centers.getAt(index); 
               iclusterIndexes = (obj.centerIndexes == index);
               subPointIndexes = obj.pointIndexes(iclusterIndexes);
               subDistances = obj.distances(subPointIndexes);
               subCenterIndexes = ones(size(subPointIndexes));
               subClusterVec = ClusterVector();
               subClusterVec.initWithFunctionSet(obj.F,subPointIndexes,subDistances,subCenter,subCenterIndexes);
        end
        
        % indexes -- local indexes
        function [subClusterVec] = subset(obj,indexes)
            subPointIndexes = obj.pointIndexes(indexes);
            % use getAt only if centers is not array 
               subDistances = obj.distances(indexes);
               subCenterIndexes = obj.centerIndexes(indexes);
               subClusterVec = ClusterVector();
                   subClusterVec.initWithFunctionSet(obj.F,subPointIndexes,subDistances,obj.centers,subCenterIndexes);
        end
        
        % Clean the object. Release memory.
        function clear(obj)
            clear obj.distances;
            clear obj.pointIndexes;
            clear obj.centerIndexes;
            clear obj.centers;
            clear obj.F;
        end
        
        function result = get.nPoints(obj)
            result = length(obj.pointIndexes);
        end
        function result = get.nClusters(obj)
            result = obj.centers.size;
        end
        
        function result=getCenterIndexes(obj, pointIndexes)
            if nargin>1
                result=obj.centerIndexes(pointIndexes);                
            else
                result=obj.centerIndexes(1:size(obj.pointIndexes,2));                
            end
            
            % ensure that center indexes is a column vector
            if size(result,1) == 1
              % should not run if we fixed it at initialization
              %mydebug
              result = result';
            end
            
        end % get the centerIndexes vector

        function result = get.d(obj)
            result = obj.F.dim;
        end
        function result=getCost(obj, costMethod)
            if costMethod == obj.maxDistanceCost 
                result=max(obj.distances);
            elseif costMethod == obj.sumDistanceCost
                result = sum(obj.distances);
            elseif costMethod == obj.sumSqDistanceCost
                result = sum(obj.distances.^2);
            else
                error('not implemented');
            end % if
        end % function
        function result=remainPoints(obj)
            indexes=setdiff(1:obj.F.M.n,obj.pointIndexes);
            result=obj.F.subset(indexes);
        end
         function show(obj,varargin)
            [varargin colored]= Utils.readParam(varargin, 'Colored', false);
             G=obj.getFunctionSet(false);
             if colored
                G.show(varargin{:});
                for i=1:obj.centers.size
                    %[~, idx] = obj.getClusterIndexes(i);
                    [temp, idx] = obj.getClusterIndexes(i);
                    [center, clusterIndexes] = ...
                        obj.getClusterIndexes(i);
                    H=obj.subsetPoints(clusterIndexes);
                    H.show('color',Utils.color(i));
                    center.show('Color',Utils.color(i));
                end % for
             else
                G.show(varargin{:});
                obj.centers.show();
             end;
         end
    end
end
