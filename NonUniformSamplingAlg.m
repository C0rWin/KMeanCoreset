classdef NonUniformSamplingAlg < handle
  % Compute a non-uniform sample from a set of numbers.
  % Usually, these number represents distances of points from a shape.
  % Each point is chosen with probability proportional to its distance
  % from a given object, and with a weight that is inverse propotional to
  % this distance.
  % Used in the construction of some coresets.
  properties
    % Number of points to sample
    sampleSize;
  end
  methods
    % set default values for the properties of the class
    function setDefaults(obj)
      Utils.completeProps(obj, {10},{'sampleSize'});
    end % function setDefaults

    % Constructor.
    % varargin (optional) -initial value for sampleSize.
    function obj=NonUniformSamplingAlg(varargin)
      if nargin>0
        obj.sampleSize=varargin{1};
      end % if
    end % function AdaptiveSamplingAlg

    function [weights, indexes, prop]= sample(obj, dists)
      % Sample indexes with probability proportional to their dist values,
      %  and assign weights proportional to the inverse of this prob.
      % Input:
      %       dists - a vector of n numbers, usually represent distance of points from
      %               a shape (usually opt center)
      %
      % Output:
      %     indexes - a vector of size sampleSize, contains the index of the chosen
      %               points (corresponding to dists)
      %     weights - a vector of sampleSize weights, where weight[i] is the
      %               weight of the i'th point

      dists(dists < 0) = 0;
      sumDists=sum(dists);

      if sumDists < eps
        mydebug
        weights = [];
        indexes = [];
        prop = [];
        return
      end

      if (nargout < 3)
        % Represent the probabilities as consecutive sub-intervals on
        % the interval [0,1], such that the length of the i'th sub-interval
        % is proportional to the i'th probability.
        % Then sample a point uniforly at random from [0,1], and choose
        % the i'th index that corresponds to the sub-interval that contains
        % this point.
        % Repeat sampleSize times, and denote by bins[i] the number of times that the
        % i'th sub-interval was selected, for every 1<i<sampleSize.

        ind = randsample(length(dists),obj.sampleSize,true,full(dists));

        % Remove the duplicated indexes
        ind= sort(ind); % can use count sort in O(n) time.
        h = histc(ind, ind);
        indexes = ind(h > 0);
        copies=h(h>0); % how many times each index appeared

        if ~isempty(obj.sampleSize)
          weights=(sumDists/obj.sampleSize)./dists(indexes);
          weights=copies.*weights;
        else
          weights = [];
        end

        if obj.sampleSize ~= sum(copies)
          error('oops, wrong number of samples');
        end
      else
        p = rand(obj.sampleSize,1)*sumDists;
        distrib = [0;cumsum(dists)];
        [~, indexes] = histc(p, distrib);
        prop = (p - distrib(indexes))./dists(indexes);
        weights=(sumDists/obj.sampleSize)./dists(indexes);
      end


    end % function sample
  end % methods

end % class NonUniformSamplingAlg

