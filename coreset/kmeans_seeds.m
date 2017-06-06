function C = kmeans_seeds(X, k, w)
    C = X(:,1+round(rand*(size(X,2)-1)));
    L = ones(1, size(X,2));
    for i = 2:k
        D = X-C(:, L);
        bsxfun(@mtimes, D, w);
        D = cumsum(sqrt(dot(D, D, 1)));
        if D(end) == 0, C(:, i:k) = X(:, ones(1, k-i+1)); return; end
        C(:, i) = X(:,find(rand < D/D(end), 1));
        [~, L] = max(bsxfun(@times, bsxfun(@minus,C'*X,dot(C,C,1)'/2), w),[],1);
    end
end
