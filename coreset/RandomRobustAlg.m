classdef RandomRobustAlg < HandleObject
%--------------------------------------------------------------------------
%This class calculates the robust median of a given set of points (with
%weights) using uniform random sample from this set. The "compute" function is the main function which fulfils the
%function. 
%
%The basic idea of the current implementation is very simple. In every
%iteration it randomly samples \beta (given as a parameter) points as the
%centers (robust median) with probability proportional to their respective
%weights. Then for every point, it calculates the distance between the
%point and its nearest center among the \beta centers. Then it sums up all
%the distances (with outliers uncounted). If it obtains smaller value, this
%set of centers will replace the previously found set of centers.
%--------------------------------------------------------------------------
    properties 
        % median fraction. Alias to medianFraction defined in AbstractBicriteriaAlg
         gamma;

         % Type = integer
         % Testable
         beta;
         figure;
         costMethod; % how to estimate cost of the centers
         nIterations;
         onSampleAlg = []; % Algorithm for running on the random sample. If gamma=1 should support compute(F) method.
         % otherwise compute(F,gamma). Default: take random sample as
         % centers
         
         partitionFraction;  % fraction of points to approximate from original set. Usually updated by bicriteria
         
         pickEndPoints; % Always choose first and last point of input set to the sample. Used for time regression problem
    end
    
    methods
        function new=copy(obj)
         new=RandomRobustAlg();
         new.gamma=obj.gamma;
         new.beta=obj.beta;
         new.figure=obj.figure;
         new.costMethod=obj.costMethod;
         new.nIterations=obj.nIterations;
         new.onSampleAlg =obj.onSampleAlg; 
         new.partitionFraction =obj.partitionFraction;
         new.pickEndPoints = obj.pickEndPoints; 
        end
        
        % compute (beta, gamma,eps)-robust median with high probability,
        % by computing gamma approximation for a random sample of beta
        % points from the original set F
        % partitionFraction is the fraction of F that
        % should approximated by the centers of sample. Needed only if
        % nIterations>1
        function optCluster = compute (obj, F)
            optCost = inf;           
            for i=1:obj.nIterations
                G=F.weightedRandomSubset(obj.beta);
                
                if obj.pickEndPoints
                    G=G.merge(F.subset(F.M.n)); % add last input point to G
                    M=F.subset(1);
                    G=M.merge(G);% add first input point to G
                    % Note that the first/last point may appear twice
                end
                if obj.figure.sample
                    G.show('clear',true,'pause', true, 'title', 'sample','pause',true)
                end
				% Generate gamma medians for the sample G.
                % if there is no robustAlg, we assume that G is the centers
                % c is a cluster
                    if obj.gamma==1
                        c = obj.onSampleAlg.computeClusterVec(G);
                    else
                        error('not implemenmted');
%                        c = obj.robustAlg.compute(G, obj.gamma);
                    end % if g=1
                if obj.figure.sample
                    if ~isa(c, 'ClusterVector')
                        c=c.getClusterVector(G);
                    end
                    Utils.show('clear',true);
                    c.show();
                    Utils.show('title', 'Sample and its robust approx');
                    Utils.show('pause',true);
                end
                % create the cluster of closest points
                % TODO: no need to copy functions!
                origCluster = ClusterVector();
                origCluster.initWithCenters(F,c.centers);
                cluster = origCluster.smallest(obj.partitionFraction);
                cost     = cluster.getCost(obj.costMethod);
                    
                % Test how good is the approximation by this iteration 
                    % compare to optimal so far
                if (i==1 || cost<optCost)
                         optCluster=cluster;
                         iOpt=i;
                         optCost=cost;
                end
                if obj.figure.iteration
                    Utils.show('clear',true);
                    R=cluster.remainPoints();
                    R.show('color','g');
                    cluster.show();
                    Utils.show('title', ['Robust iteration ' num2str(i) ': Removed partitionFraction and their approximation']);
                    R.show('clear',2);
                    origCluster.show('Colored',true,'clear',2);
                    Utils.show('pause', true);
                end  % if figure
            end % for
              if obj.nIterations > 1
                if obj.figure.opt
                        F.show('color','g','clear',true);
                        cluster.show();
                        Utils.show('pause', true, 'title', ['Optimal Robust: ' num2str(iOpt)]);
                end % if figure
              end
        end % function compute

    end % methods
end % class 