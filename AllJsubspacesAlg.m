classdef AllJsubspacesAlg < HandleObject
    properties
        j=0;
        regression=false;
        throughOrigin=false;
    end
    
    methods
        function c=computeClusterVec(obj, F)
             c = obj.compute(F);
        end
        % compute all j-subspaces that are spanned by j points of F.matrix
        % Distances to F are not updated in clusterVec
        function c = compute(obj,F) % 
            if obj.j==0
                c = ClusterVector();
                c.initWithFunctionSet(F, 1:F.M.n, zeros(F.M.n,1), F, 1:F.M.n)

            else
                subspaces=Subspaces();
                subspaces.allJSubspaces(obj.j, F, obj.regression, obj.throughOrigin);
                c = ClusterVector();
                c.initWithCenters(F,subspaces); 
            end
        end % subspaces
    end
    
end

