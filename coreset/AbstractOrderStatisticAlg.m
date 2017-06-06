classdef AbstractOrderStatisticAlg
%--------------------------------------------------------------------------
%This abstract class is the base class for all classes which are to compute
%some statistical quantity according to some sorts. Of course a complete
%sorting is not necessary in order to figure out certain quantity.
%--------------------------------------------------------------------------
     methods (Abstract)
         % mValue: the number such that the sum of elements of "vector"
         % whose "vals" value is less than or equal to mValue is at least
         % m.
         % indexes: (NOT-NECESSARILY-SORTED) indexes of the "vals" vector
         % whose values are less than or equal to mValue. 
         [indexes, mValue] = compute (vals, m, vector, minimum);
     end
end

