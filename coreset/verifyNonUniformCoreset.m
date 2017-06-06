function verifyNonUniformCoreset( )

    N = 50000;
    
    D = 10;

    P = rand(N, D);

    P = bsxfun(@minus, P, mean(P, 1));
    
    W = ones(1, N);
    
    k = 1;
    
    [C, L] = kmeanspp_fast_weighted(P', k, W);
    
    distance = zeros(1, N);
    for c=1:k
        distance(:,(L==c)') = sum(bsxfun(@times, bsxfun(@minus, P((L==c)',:), C(:,c)').^2, W(L==c)'), 2);
    end
    distance = distance';

    % Compute sensitivity function based on approximated clusters for each
    % point
    s_p = zeros(N, 1);
    for c=1:k        
        s_p((L==c),:) = W(L==c)'*8./sum(W(L==c),2) + 2*distance(L==c).*W(:, L==c)'./sum(distance);
    end
    
    [~, optCenters] = Ckmeans(P, k, ones(N, 1), 'distance', 'sqeuclidean', 'maxiter', 1000);
    optEnergy = computeEnergy(P, optCenters, k);
    
    i = 0;
    % Running tests where parameter t
    % defines the coreset size, trying to build coresets
    le = length(1000:100:20000);
    error = zeros(le, 1);
    
    for t=1000:100:20000
        i = i + 1;

        % compute points weights based on sensitivity function and sample
        % size. w(i) is the weight of ith point in dataset
        weights = bsxfun(@rdivide, (W*(sum(s_p)))', (s_p*t));

        % sample t points
        sample_idx = randsample(N, t, true, s_p);

        % get sampled points (desired coreset)
        coreset = P(sample_idx,:);

        % calculate coreset center
        [~, ccenter] = Ckmeans(coreset, k, weights(sample_idx), 'distance', 'sqeuclidean', 'maxiter', 1000);
        % coresetEnergy = computeEnergy(P, ccenter, k);
        
        % calculate error, based on weak coreset definition
        % abs(coresetEnery - optEnergy) <= error * optEnergy
        % cost(P,opt(C))~cost(P,OPT(P))
        
        %error(i) = coresetEnergy/optEnergy - 1;
        error(i) = norm(ccenter);
    end

    figure
    scatter(1:length(error), error, 'b*');
    xlabel('size of T');
    ylabel('error');
    legend('Random sampling');
end

% cost function
function energy = computeEnergy(points, centers, k)
    dist = zeros(size(points,1), k);
    for c=1:k
        dist(:,c) = sum( bsxfun(@minus, points, centers(c,:)).^2, 2);
    end
    energy = sum(min(dist, [], 2));
end
