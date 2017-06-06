function [ centers, labels ] = my_plusplus( X, k )

    [n, d] = size(X);

    labels  = ones(n, 1);

    centers = zeros(k, d);

    centers(1, :) = X( randi(n),:);

%     while length(unique(labels) ~= k)
        for i=2:k
            dist = X - centers(labels, :);
            dist = dot(dist, dist, 2);
            prob = dist./sum(dist);
            centers(i,:) = X(randsample(1:n, 1, true, prob), :);
            [~, labels] = max(bsxfun(@minus, X*centers', dot(centers, centers, 2).'/2), [], 2);
        end
%     end
end
