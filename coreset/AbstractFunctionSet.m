classdef AbstractFunctionSet  < HandleObject
%--------------------------------------------------------------------------
%This class is an abstract class for a set of (weighted) points or a
%subspace.
%--------------------------------------------------------------------------
    properties
      % Type=Matrix
      % A matrix of size n-by-d. Used to store a set P of n points in R^d
      M;

    end 
    
    properties (Dependent)
      
      % Type = n
      size;
      
      %type = integer
      dim;
    end
        
    methods (Abstract)
        % obj is AbstractFunctionSet. Need downcast/template.
        vals = Eval(obj, center); %vals = Eval(AbstractFunctionSet obj, Subspace center)
        G = subset(obj, indexes); %Return the subset of obj consists of row indexed as indexes.
        S = randomSubset(obj, sampleSize, repetition); %Return a random subset of obj with size sampleSize
    end
   
    methods
      % Constructor
      function obj = AbstractFuntionSet(FunctionSet)
          obj.M = FunctionSet;
      end
      
      function result = get.dim(obj)
          result = obj.M.nCols;
      end 
      
      function result = get.size(obj)
          result = obj.M.nRows;
      end
    end
end 