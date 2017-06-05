classdef PointFunctionSet < AbstractFunctionSet
%--------------------------------------------------------------------------
%This class represents a set of weighted points where each point is a
%function.
%--------------------------------------------------------------------------

    properties
        % Weights of the points. Vector of n elements.
        % NOTE: W is n-by-1 vector.
        % Type = Matrix
        W=Matrix(); 
        
        % Type =  real
        com; % center of mass for addDim
        highVal=realmax;
    end
    properties (Dependent)
        sumW;
        n;
    end % properties (Dependent)
    
    methods
        function [Q idxs]=epsilonGrid (obj, epsilon, G)
            %[~, idxs] = obj.M.epsilonGrid(epsilon, G.M);
            [temp, idxs] = obj.M.epsilonGrid(epsilon, G.M);
            Q=G.subset(idxs);
        end
        function out=getAt(obj, i)
            out=obj.subset(i);
        end
        function result=get.n(obj)
            result = obj.M.n;
        end
        function result=get.sumW(obj)
            if obj.W.isEmpty()
                result=obj.M.n;
            else
                result=sum(obj.W.m);
            end
        end
        
        % Instantiate new object of the same class.
        function new=copy(this)
           new = feval(class(this));
           new.copyToExistent(this);
        end
        function copyToExistent(obj, other)
           obj.M=Matrix(other.M.m);
           obj.W=Matrix(other.W.m); 
        end
        function addDim(obj)
            m=obj.M.m;
            obj.com=mean(m);
            newM=Matrix(bsxfun(@minus, m, obj.com));
            obj.highVal = newM.sumOfSquaredEntries()*10^8;
            obj.M = Matrix([newM.m obj.highVal*ones(obj.M.n,1)]);
%            obj.M = Matrix(bsxfun(@plus, obj.M.m(:,1:end-1), obj.com));
        end
        
        function F=removeDim(obj)
            m=obj.M.m;
            last=obj.highVal./m(:,end);
            newM=bsxfun(@times,m(:,1:end-1),last);
            obj.M = Matrix(bsxfun(@plus, newM, obj.com));
            F=obj;
        end
        
        function result=isEmpty(obj)
            result= obj.M.n==0;
        end
        % Constructor
        % Points is a Matrix object (or empty set)
        function obj = PointFunctionSet(Points, varargin)
            if nargin == 0
                obj.M = Matrix();
                obj.W = Matrix();
            end
            if nargin == 1
                if isa(Points, 'double')
                    Points=Matrix(Points);
                end
                obj.M = Points;
                obj.W = Matrix(ones(obj.size, 1));
            end
            if nargin > 1
                obj.M = Points;
                obj.W = varargin{1};
            end
        end
               
        % Compute squared distance from obj (set of points) to centers
        % Either set of points or cell of subspaces        
		% Both obj and centers are PointFunctionSet type.
        % We also return the indices that reach the minimum.
        
        function [m argm]=maxDistance(obj,centers)
            sqDistances = obj.Eval(centers);
            [m2 argm]=max(sqDistances);
            m=sqrt(m2);
        end            
        
        function m=sumDistance(obj,centers)
            sqDistances = obj.Eval(centers);
            m=sum(sqrt(sqDistances));
        end  
        % This function should just call centers.Eval.
        function [sqDistances indexes c] = Eval(obj, centers)
            if isa(centers, 'GISPoints')
                F=obj.withoutTime;
                centers=centers.withoutTime;
            else
                F=obj;
            end
            if isa(centers,'PointFunctionSet')
                N=NearestCenterAlg();
                [sqDistances indexes] = N.compute(F.M.matrix, F.W.matrix, centers.M.matrix);
            elseif isa(centers,'Subspace')
                % Compute squared distance from point in Matrix to subspace
                sqDistances= obj.M.sqDistances(centers);
            
                %indexes is created to match for the abstract Eval, in clusters
                %there are indexes for the nearest pointFunctionSet
                indexes = ones(size(sqDistances,1),size(sqDistances,2));
            elseif isa(centers,'Subspaces')
                % Compute squared distance from point in Matrix to subspace
                [sqDistances indexes] = centers.sqDistances(obj);
                c=ClusterVector();
                c.initWithFunctionSet(obj, 1:obj.size, sqDistances, centers, indexes);
            end
        end
        
        function setSubset(obj, indexes, G)
            try
                obj.M.m(indexes,:)=G.M.m;
                obj.W.m(indexes,:)=G.W.m; 
            catch
                error('df')
            end
        end
        
        % Return the subset (new copy) of matrix M consists of rows indexed by indexes
        function G = subset(obj, indexes)
          G=obj.copy();
          G.M=obj.M.getRows(indexes);
          G.W=obj.W.getRows(indexes);
        end
        
        % Return a random subset of obj with size sampleSize
        function  S = randomSubset(obj, sampleSize, repetition)
          [randM, randIndex] = obj.M.sampleRows(sampleSize,repetition);  
          S = PointFunctionSet(randM, obj.W.getRows(randIndex));
        end
        
        %  Return a
        % subset of points in G with size sampleSize as the medians. The points are
        % sampled according to their weights in G.
        % copies says how much time we selected the corresponding point in
        % sample
        function [sample copies]=weightedRandomSubset(obj, sampleSize)
            % Initialize a non-uniform sampler
            nonUniformSampler = NonUniformSamplingAlg(sampleSize);
            [weights,indexes] = nonUniformSampler.sample(obj.W.matrix);
            copies = weights*sampleSize/sum(weights);
            sample=obj.subset(indexes);
        end
            
        % Remove certain subset of points according to the indexes
        function result = remove(obj, mIndexes)
            rest = 1:obj.size;
            rest(mIndexes)=[];
            result = obj.subset(rest);
        end
        
        % Merge two PointFunctionSets
        function obj = merge(obj, P2)
            obj.M.merge(P2.M);
            obj.W.merge(P2.W);
%             % can not use merge function from Matrix since we need to check
%             % dimensions
%             if obj.M.nCols==0
%                 % obj.dim = P2.dim;
%                 obj.M.matrix=[obj.M.matrix; P2.M.matrix];
%                 obj.W.matrix=[obj.W.matrix; P2.W.matrix];
%                 obj.sumW = sum(obj.W.matrix);
%                 % obj.size = obj.M.nRows;
%             elseif P2.M.nCols== 0
%                 % Do nothing
%             elseif obj.M.nCols == P2.M.nCols
%                 %
%                 obj.M.matrix=[obj.M.matrix; P2.M.matrix];
%                 obj.W.matrix=[obj.W.matrix; P2.W.matrix];
%                 obj.sumW = sum(obj.W.matrix);
%                 % obj.size = obj.M.nRows;
%             else
%                 warning('Cannot merge pointsets with different in dimensions');
%             end
        end %end merge
        
        function clear(obj)
           obj.M.clear();
           obj.W.clear();
           obj.sumW = 0;
        end
      
        function subspace = getSpannedSubspace(obj)
            subspace = obj.M.getSpannedSubspace();
        end
        function result = getRows(obj,indexes)
            result=obj.copy();
            result.M=obj.M.getRows(indexes);
            result.W=obj.W.getRows(indexes);
        end


        function show(obj, varargin)
            Utils.show('P', obj, 'func', 'plot3' ,'.',varargin{:});
        end % show
                
    end % methods
end