function [centers, label] = kmeanspp_fast_weighted(X, k, w)
m = seeds(X, k, w);
%[A, U] = bicriteria(X', k, w', 55);
%m = seeds(A', k, U');
[label, centers] = kmeans(X, m, w);
%idx = knnsearch(X', y', 'k', 1);
%centers = X(:, idx);

function m = seeds(X, k, w)
m = X(:,1+round(rand*(size(X,2)-1)));
L = ones(1, size(X,2));
for i = 2:k
    D = X-m(:, L);
    bsxfun(@mtimes, D, w);
    D = cumsum(sqrt(dot(D, D, 1)));
    if D(end) == 0, m(:, i:k) = X(:, ones(1, k-i+1)); return; end
    m(:, i) = X(:,find(rand < D/D(end), 1));
    [~, L] = max(bsxfun(@times, bsxfun(@minus,m'*X,dot(m,m,1)'/2), w),[],1);
end

function [label, centers] = kmeans(X, m, weights)
n = size(X,2);
last = 0;
[~,label] = max(bsxfun(@times, bsxfun(@minus,m'*X,dot(m,m,1)'/2), weights),[],1);
while any(label ~= last)
%for i=1:5
    [u,~,label(:)] = unique(label);
    k = length(u);
    E = sparse(1:n, label, weights, n, k, n);
    m = X*(E*spdiags(1./sum(E,1)', 0, k, k));
    last = label;
    [~, label] = max(bsxfun(@times, bsxfun(@minus,m'*X,dot(m,m,1)'/2), weights),[],1);
end
centers = m;