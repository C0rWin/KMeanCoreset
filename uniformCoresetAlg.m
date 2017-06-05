classdef uniformCoresetAlg

    properties
        t;
    end

    methods
        % Merge two coresets C1, C2 to a new one.
        % C1 and C2 are PointfunctionSets. Compute the coreset on the merged
        % set of points. This is a merged Coreset.
        function C=mergedCoreset(obj, C1, C2)
          % Coresets are PointFunctionSets
          C1.merge(C2);
          C = obj.computeCoreset(C1);
        end

        function c = computeCoreset(obj, P)
            try
                coreset = datasample(obj.matrix, obj.t);
                weights = ones( obj.t,1)*P.size/ obj.t;
                c = matrixToFunctionSet(coreset, weights);
            catch err
                disp(err);
            end
        end

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
