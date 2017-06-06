    classdef nonUniformCoreset
    properties
        k;
        t;
        maxiter;
    end

    methods

        function obj = nonUniformCoreset(k, t, maxiter)
            obj.k = k;
            obj.t = t;
            obj.maxiter = maxiter;
        end

        function [s_p] = sensitivityWithSeed(obj, P, L, C)
            distance = zeros(1, P.size);
            for c=1:obj.k
                distance(:,(L==c)') = sum(bsxfun(@minus, P.M.m((L==c)',:), C(:,c)').^2,2);
            end
            distance = distance';
            s_p = zeros(P.size, 1);
            for c=1:obj.k
                s_p((L==c),:) = 1/size(L==c,2) + distance((L==c)')./sum(distance);
            end
        end

        function [s_p] = sensitivity(obj, P)
            [L, C] = kmeanspp_weighted(P.M.m', obj.k, P.W.m);
            distance = zeros(1, P.size);
            for c=1:obj.k
                distance(:,(L==c)') = sum(bsxfun(@minus, P.M.m((L==c)',:), C(:,c)').^2,2);
            end
            distance = distance';
            s_p = zeros(P.size, 1);
            for c=1:obj.k
                s_p((L==c),:) = 1/size(L==c,2) + distance((L==c)')./sum(distance);
            end
        end

        function C = mergedCoreset(obj, C1, C2)
          % Coresets are PointFunctionSets
          C1.merge(C2);
          C = obj.computeCoreset(C1);
        end

        function c = computeCoreset(obj, P)
            try
                %take care of weighted input
                weights = P.W.m./(obj.sensitivity(P)*obj.t);
                if P.size < obj.t
                    % c = obj.matrixToFunctionSet(P.M.m, weights);
                    c = obj.matrixToFunctionSet(P.M.m, ones(P.size, 1));
                else
                    sample_idx = randsample(P.size, obj.t, true, obj.sensitivity(P));
                    coreset = P.M.m(sample_idx,:);
                    c = obj.matrixToFunctionSet(coreset, weights(sample_idx));
                end;
            catch err
                disp(err);
            end
        end

    end

    methods (Static)
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
