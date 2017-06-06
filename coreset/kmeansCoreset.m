classdef kmeansCoreset

    properties
        t;
        maxiter;
    end

    methods

        function obj = kmeansCoreset(t, maxiter)
            obj.t = t;
            obj.maxiter = maxiter;
        end

        function C = mergedCoreset(obj, C1, C2)
          % Coresets are PointFunctionSets
          try
            C1.merge(C2);
%             if (size(C1.M.m,1) ~= size(C1.W.m,1))
%                 disp('Err!');
%             end
            C = obj.computeCoreset(C1);
          catch err
              err.getReport;
          end
        end

        function c = computeCoreset(obj, P)
            % kmeanspp_fast_weighted(M', i*factor, ones(1, size(M, 1)))
            try
                if P.size < obj.t
                    c = obj.matrixToFunctionSet(P.M.m, ones(P.size, 1));
                else

                     [t_centers, m_part] = kmeanspp_fast_weighted(P.M.m', obj.t, P.W.m');
                     t_centers = t_centers';

                     %[m_part, t_centers, ~, ~] = Ckmeans(P.M.m, obj.t, P.W.m,'distance', 'sqeuclidean', ...
                     %   'maxiter', 100, 'emptyaction', 'singleton', 'display', 'off', 'onlinephase', 'off');

                    weights = zeros(obj.t, 1);
                    for i=1:obj.t
                        weights(i) = sum(P.W.m(m_part==i));
                    end
                    c = obj.matrixToFunctionSet(t_centers, weights);
                end
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
