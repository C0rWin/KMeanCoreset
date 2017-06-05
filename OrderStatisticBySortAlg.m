classdef OrderStatisticBySortAlg < AbstractOrderStatisticAlg
%--------------------------------------------------------------------------
%This class is to figure out some statistical quantity according to some
%sorts. Of course in the implementation we don't really need to do the
%sorting itself.
%--------------------------------------------------------------------------
    methods (Static)
        % Input: numbers - vector of number to find the order statistics from
        %        os - the order statistics (between 0 to sum of weights)
        %        weights - weight vector (the same size of 'numbers'). Weights
        %        are positive reals.
        % 
        
        % Output: Minimum number of points to return
        % mValue - the weighted median number of 'numbers'. This is the number
        %          where the sum of weights of points with numbers' value less
        %          than or equal to mValue is at least os.          
        % indexes - (NOT-NECESSARILY-SORTED) indexes of the numbers with value less than or
        %           equal to mValue.

        function [indexes mValue]=compute(numbers, os, weights, minimum)
           %[~,I]=sort(numbers);
           [B,I]=sort(numbers);
           if nargin < 4
               minimum = 1;
           end
           if nargin == 2
               m=os;
           else
               % Figure out the index m in I such that the sum of weights of
               % the points in indexes I(1:m) is os
               cumW = cumsum(weights(I));
               m = find(cumW >= os, 1, 'first');
               if m < minimum
                   m = minimum;
               end
           end
           try
               indexes=I(1:m); %os smallest numbers' indexes
           catch
               warning('debug: OrderStatisticAlg got minimum index larger then set size');
               indexes=I;
               m=length(I);
           end
           mValue = numbers(I(m)); %the os value
        end 
    end   
end
