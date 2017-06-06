function energy = computeEnergy(points, centers, k)
    dist = zeros(size(points,1), k);
    for c=1:k
        dist(:,c) = sum( bsxfun(@minus, points, centers(c,:)).^2, 2);
    end         
    energy = sum(min(dist, [], 2));
end