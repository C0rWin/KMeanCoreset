classdef uniformCoreset
    properties
        t;
    end

    methods

        function obj = uniformCoreset(t)
            obj.t = t;
        end

        function C = mergedCoreset(obj, C1, C2)
          % Coresets are PointFunctionSets
          C1.merge(C2);
          C = obj.computeCoreset(C1);
        end

        function c = computeCoreset(obj, P)
            try
                if P.size < obj.t
                    c = obj.matrixToFunctionSet(P.M.m, ones(P.size, 1));
                else
                    coreset = datasample(P.M.m, obj.t, 'Weights', P.W.m);
                    c = obj.matrixToFunctionSet(coreset, ones(obj.t, 1)*P.size/obj.t);
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
