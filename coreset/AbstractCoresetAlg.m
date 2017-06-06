classdef AbstractCoresetAlg < handle

    properties
    % Type=AbstractBicriteriaAlg 
    % Used to compute bi-criteria and then apply the coreset on each 
    % cluster of the output, using computeUsingBicriteria().
    bicriteriaAlg;
    end
    
    % Abstract class that implements functions needed for testing
    % the coreset (QA), and for use with the streaming class.
    %   Detailed explanation goes here
    methods
        function setDefaults(obj)
        end
        
        %compute the coreset either with bicriteria or regular
        %input: P - Type=AbstractFunctionSet, Subspace
        function coreset = computeCoreset(varargin) %(obj,P,subspace)
            obj = varargin{1};
            if isempty(obj.bicriteriaAlg)
                coreset = obj.compute(varargin{2:end});
            else
                coreset = obj.computeUsingBicriteria(varargin{2});
            end
        end
        
        % Compute bicriteria approximtion and apply the coreset on each
        % cluster. Returns the union of coresets.
        function unitedCoreset=computeUsingBicriteria(obj,F, bicriteriaAlg)
            if isempty(obj.bicriteriaAlg)
                obj.bicriteriaAlg = bicriteriaAlg;
            end
            clusters=obj.bicriteriaAlg.compute(F,1);
            for nCluster=1:clusters.centerCounter
                currCluster=clusters.getCluster(nCluster); %need a function for subset
                coreset=obj.computeUsingCluster(currCluster);
                if nCluster==1
                    unitedCoreset=coreset;
                else
                    %merging all the coresets in one operation will
                    %optimise the running time
                    unitedCoreset.merge(coreset);
                end % if
            end % for
        end % computeUsingBicriteria

        % Merge two coresets C1, C2 to a new one.
        % Output a single coreset of the union of C1 and C2.
        function C=mergedCoreset(obj, C1, C2)
            C1.merge(C2);
            if isempty(obj.bicriteriaAlg)
                C = obj.computeCoreset(C1);
            else
                C=obj.computeUsingBicriteria(C1,obj.bicriteriaAlg);
            end
        end % function

        % Sanity QA check for a coreset.
        % Compute a random matrix P using varargin, and check that the
        % optimum of P is approximately the optimum of its coreset C
        % varargin - see constructor of Matrix.
        function  [errorCostPUsingOptC, errorOptCostC ] = testRegular(obj, varargin)
            % if no argument used, at least set n
            obj.setDefaults();
            if nargin==1
                varargin ={'dense'};
            end % if
            P =Matrix(varargin{:});
            C=computeCoreset(obj, P);
            [errorCostPUsingOptC, errorOptCostC] =  ...
             obj.isGoodCoreset(P, C);
        end
        
       
        % Construct two independent random matrices P1,P2 by calling
        % Matrix(varargin).
        % Then, compute corresponding coresets C1,C2, and check whether the
        % optimum of the merged coreset is the same as the optimum of
        % P1, P2.
        % Such a merged step is used by the stream class.
        function  [errorCostPUsingOptC, errorOptCostC] = testMerge(obj, varargin)
            obj.setDefaults();
            if nargin==1
                varargin={'dense'};
            end
            P1 =Matrix(varargin{:});
            P2 =Matrix(varargin{:});
            C1=obj.computeCoreset(P1);
            C2=obj.computeCoreset(P2);
            C=mergedCoreset(obj,C1, C2);
            [errorCostPUsingOptC, errorOptCostC]=obj.isGoodCoreset(P1.merge(P2), C);
        end
        
        % Sanity QA check for a coreset.
        % Compute a random matrix P using varargin, and check that the
        % optimum of P is approximately the optimum of its coreset C
        % varargin - see constructor of Matrix.
        function [errorCostPUsingOptC, errorOptCostC]=isGoodCoreset(obj, P, C)
            % Check that C is a weak coreset for P.

            % Compute the optimum  of P
            optCostForP=obj.computeOptCost(P);
            
            % Compute the optimum of C.
            optCostForC=obj.computeOptCost(C);
            
            % Check that they have approximately the same cost.
            errorOptCostC=obj.relativeError(optCostForP,optCostForC);

            % Compute the cost from P to the optimum of the coreset,
            % compare to the real opt of P. 
            % Unlike the previous check, here the optimum of the coreset
            % is applied on P, not on C.
            optForC=obj.computeOpt(C);
            costPUsingOptC=obj.computeCost(P, optForC);
            errorCostPUsingOptC=obj.relativeError(optCostForP,costPUsingOptC);
        end % function isGoodCoreset
        function [C time]=computeCoresetWithTime(obj, P)
            tic;
            C=obj.computeCoreset(P);
            time=toc;
        end
        
    end % methods
%     methods (Abstract)        
%        % Compute the optimal cost of the set P
%        optCost=computeOptCost(obj, P);
% 
%        % copmute the cost of the center Q to the set P
%        cost=computeCost(obj, P, Q);
% 
%        % Compute the optimal solution of P
%        opt=computeOpt(obj, P);
% 
%        % Compute a coreset for P
% %        C=compute(obj, P);
%     end % abstract methods
    
    methods (Access=protected, Static)
        function error = relativeError (opt, approx)
            % return a number between 0 to 1, represents the relative error of
            % approximating approx by opt
            error=Utils.relativeError(opt,approx);
        end % function error=..
    end %  static methods
end % class AbstractCoresetAlg


