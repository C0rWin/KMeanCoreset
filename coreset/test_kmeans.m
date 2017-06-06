function test_kmeans()

    warning off;

    load data/mnist_all

    M = double([train0; train1; train2; train3; train4; train5]);    
    
    trials = 1:5;

    factor = 5;
    
    indexes = factor * trials;
    
    T = zeros(7, length(trials));
    
    E = zeros(7, length(trials));

    for i=trials, tic; [c, ~] = kmeanspp_fast_weighted(M', i*factor, ones(1, size(M, 1))); T(1, i) = toc; E(1, i) = computeEnergy(M, c', size(c, 2)); end;       
    
    for i=trials, tic; [~, c, ~, ~] = kmeans(M, i*factor, 'MaxIter', 1); T(2,i) = toc; E(2, i) = computeEnergy(M, c, i*factor); end;
    
    for i=trials, tic; [A, U] = bicriteria(M, i*factor, ones(size(M, 1), 1), 60); [c, ~] = kmeanspp_fast_weighted(A', i*factor, U'); T(3, i) = toc; E(3, i) = computeEnergy(M, c', size(c, 2)); end;
    
    for i=trials, tic; [A, U] = bicriteria(M, i*factor, ones(size(M, 1), 1), 75); [c, ~] = kmeanspp_fast_weighted(A', i*factor, U'); T(4, i) = toc; E(4, i) = computeEnergy(M, c', size(c, 2)); end;
    
    for i=trials, tic; [A, U] = bicriteria_and_seed(M, i*factor, ones(size(M, 1), 1), 60); [c, ~] = kmeanspp_fast_weighted(A', i*factor, U'); T(5, i) = toc; E(5, i) = computeEnergy(M, c', size(c, 2)); end;       
    
    for i=trials, tic; [A, U] = bicriteria_and_seed(M, i*factor, ones(size(M, 1), 1), 75); [c, ~] = kmeanspp_fast_weighted(A', i*factor, U'); T(6, i) = toc; E(6, i) = computeEnergy(M, c', size(c, 2)); end;       
    
    for i=trials, tic; c = kmeans_seeds(M', i*factor, ones(1, size(M, 1))); T(7, i) = toc; E(7, i) = computeEnergy(M, c', size(c, 2)); end;       
    
    % Error
    figure;
    plot(indexes, E(1,:), 'b-o', indexes, E(2,:), 'r:o', ...
        indexes, E(3, :), 'g-.d', indexes, E(4, :), 'm-.x', ...
        indexes, E(5, :), 'k-*', indexes, E(6, :), 'y-*', ...
        indexes, E(7, :), 'c-.s', 'LineWidth', 1.5);
    
    legend('kmeans++', 'matlab kmeans', 'bicriteria p=60 & kmeans++', ...
        'bicriteria p=75 & kmeans++', 'bicriteria+seed p=60 & kmeans++', ...
        'bicriteria+seed p=75 & kmeans++', 'seed');
     xlabel('Number of clusters');
     ylabel('Cost (sum of squared distances)');
     title('Cost function graph vs cluster size');
 
     % Time
     figure
     plot(indexes, T(1,:), 'b-o', indexes, T(2,:), 'r:o', ...
         indexes, T(3, :), 'g-.d', indexes, T(4, :), 'm-.x', ...
         indexes, T(5, :), 'k-*', indexes, T(6, :), 'y-*', ...
         indexes, T(7, :), 'c-.s', 'LineWidth', 1.5);
     
    legend('kmeans++', 'matlab kmeans', 'bicriteria p=60 & kmeans++', ...
        'bicriteria p=75 & kmeans++', 'bicriteria+seed p=60 & kmeans++', ...
        'bicriteria+seed p=75 & kmeans++', 'seed');
     xlabel('Number of clusters');
     ylabel('Execution time (ms)');
     title('Running times vs cluster size');

%     clear;
% 
%     warning off;
% 
%     load data/mnist_all
% 
%     M = double([train0; train1; train2; train3; train4; train5]);   
%     
%     k = 100;
%     T = [];
%     E = [];
%     
%     step = 1000;
%     trials = 1:20;
%     
%     [c_kpp, ~] = kmeanspp_fast_weighted(M', k, ones(1, size(M, 1))); 
%     E(1, :) = computeEnergy(M, c_kpp', size(c_kpp, 2)) * ones(length(trials), 1);
%     
%     [~, c_k, ~, ~] = kmeans(M, k, 'MaxIter', 1); 
%     E(2, :) = computeEnergy(M, c_k, k) * ones(length(trials), 1);
%     
%     for i=trials, tic; [A, U] = bicriteria(M, k+step*i, ones(size(M, 1), 1), 50); [c, ~] = kmeanspp_fast_weighted(A', k, U'); T(1, i) = toc; E(3, i) = computeEnergy(M, c', size(c, 2)); end;
%     
%     for i=trials, tic; [A, U] = bicriteria(M, k+step*i, ones(size(M, 1), 1), 55); [c, ~] = kmeanspp_fast_weighted(A', k, U'); T(2, i) = toc; E(4, i) = computeEnergy(M, c', size(c, 2)); end;
%     
%     for i=trials, tic; [A, U] = bicriteria(M, k+step*i, ones(size(M, 1), 1), 60); [c, ~] = kmeanspp_fast_weighted(A', k, U'); T(3, i) = toc; E(5, i) = computeEnergy(M, c', size(c, 2)); end;       
%     
%     for i=trials, tic; [A, U] = bicriteria(M, k+step*i, ones(size(M, 1), 1), 65); [c, ~] = kmeanspp_fast_weighted(A', k, U'); T(4, i) = toc; E(6, i) = computeEnergy(M, c', size(c, 2)); end;       
%     
%     for i=trials, tic; [A, U] = bicriteria(M, k+step*i, ones(size(M, 1), 1), 75); [c, ~] = kmeanspp_fast_weighted(A', k, U'); T(5, i) = toc; E(7, i) = computeEnergy(M, c', size(c, 2)); end;       
%     
%      % Error
%     figure;
%     plot(trials, E(1,:), 'b-o', trials, E(2,:), 'r:o', ...
%         trials, E(3, :), 'g-.d', trials, E(4, :), 'm-.x', ...
%         trials, E(5, :), 'k-*', trials, E(6, :), 'y-*', ...
%         trials, E(7, :), 'c--s', 'LineWidth', 1.5);
%     
%     legend('kmeans++', 'matlab kmeans', 'bicriteria p=50 & kmeans++', ...
%         'bicriteria p=55 & kmeans++', 'bicriteria p=60 & kmeans++', ...
%         'bicriteria p=65 & kmeans++', 'bicriteria p=75 & kmeans++');
%     xlabel('Bicriteria sample size');
%     ylabel('Cost (sum of squared distances)');
%     title('Cost function graph vs cluster size');
% 
%     % Time
%     figure
%     plot(trials, T(1, :), 'g-.d', trials, T(2, :), 'm-.x', ...
%         trials, T(3, :), 'k-*', trials, T(4, :), 'y-*', ...
%         trials, T(5, :), 'c--s', 'LineWidth', 1.5);
%     
%     legend('bicriteria p=50 & kmeans++', ...
%         'bicriteria p=55 & kmeans++', 'bicriteria p=60 & kmeans++', ...
%         'bicriteria p=65 & kmeans++', 'bicriteria p=75 & kmeans++');
%     xlabel('Bicriteria sample size');
%     ylabel('Execution time (ms)');
%     title('Running times vs cluster size');   
    
end
