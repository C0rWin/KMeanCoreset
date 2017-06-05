classdef  NearestCenterAlg 
    % Given a set F of n points with dimension d and a set C of m points
    % with dimension d as centers, for every point p in F, calculate its
    % squared Euclidean distance to the nearest point in C. Return the distance 
    % and the index of its nearest center in C.
    % In the weighted case, every point in F has a weight. The final
    % distance for every point would be the geometric distance times their
    % respective weights.
    
    methods 
        function [sqDistances indexes] = compute(obj, F, FW, C, method)
            if ~exist('method','var') || isempty(method)
                method='mtimes';
            end % if
            switch (lower(method))
                case 'cellfun'
                    [sqDistances indexes] = obj.computeCellFun(F, FW, C);
                case 'mtimes'
                    [sqDistances indexes] = obj.computeMtimes(F, FW, C);
                 otherwise
                    error('unkown input parameter');
            end % switch
        end % compute

        function result=myMtimes(obj, A, B)
            mulTime=tic;
            result=A*B;
            mulTime=toc(mulTime);
            %disp(['Multiplication Time: ' num2str(mulTime)]);
        end % myMtimes
        
        function [sqDistances indexes] = computeMtimes(obj, F, FW, C)
            
            % min_c |c-p|^2=min_c c*c+p*p-2*c*p =min_c c*c-2*c*p
            % P1 is necessary only for SqDistances. Probably can be
            % removed.
            C1=2*C;
            FC=obj.myMtimes(F,C1'); % entry ij is 2*p*c for the ith point and jth center

            C2=sum(C.^2,2);
            diff=bsxfun(@minus, C2', FC); % entry ij is c*c-2*p*c for the ith point and jth center
            [mdif, indexes]=min(diff,[],2); %=min(diff');
            F1=sum(F.^2,2);
            nonWeightedSqDistances=F1+mdif;
            sqDistances=FW.*nonWeightedSqDistances;
        end % computeMtimes(F, FW, C)
        
        % Do the computation as described above
        function [sqDistances indexes] = computeCellFun(obj, F, FW, C)
%       Solution 1 (Learn from MATLAB kmeans):
%             [n,d] = size(F);
%             D = zeros(n,size(C,1));
%             nclusts = size(C,1);
% 
%            for i = 1:nclusts
%                 D(:,i) = FW(:,1).* (F(:,1) - C(i,1)).^2 ;
%                 for j = 2:d
%                     D(:,i) = D(:,i) + FW(:,1).*(F(:,j) - C(i,j)).^2;
%                 end
%            end

%       Solution 2:
            center = Matrix(C);
            fun = @(c)(NearestCenter.sqDistsMatrixPoint(F, c, FW));
            centerCell = mat2cell(C, ones(center.n, 1), center.d);
            D = cellfun(fun, centerCell, 'UniformOutput', false);
            D = cell2mat(D');
            [sqDistances indexes] = min(D,[],2); 
        end
    end
    
    methods (Access = private, Static)
        function sqDists = sqDistsMatrixPoint(P, c, weight)
            distVectors = bsxfun(@minus, P, c);
            sqDists=Matrix(distVectors).squaredRows().* weight;
        end
    end
end
