classdef KmeansCoresetTest < Test

    properties

        caseNo; n; d; k; t; realSize; matrix;

        optEnergy; uniformEnergy; nonUniformEnergy; kmeansEnergy;

        uniformError; nonUniformError; kmeansError;

        kmeansErrorWithSearch;

        maxiter = 100; display = 'off';

        synthetic = true;

        fullKmeanTime, uniformTime, nonUniformTime, kmeanCoresetTime;

        isFastKmeans = false;

    end

    methods

        function energy = computeEnergy(obj, centers)
            dist = zeros(size(obj.matrix, 1), obj.k);
            for c=1:obj.k
                dist(:,c) = sum(bsxfun(@minus, obj.matrix, centers(c,:)).^2, 2);
            end
            energy = sum(min(dist, [], 2));
        end

        function [s_p] = sensitivity(obj)
            [L, C] = kmeanspp(obj.matrix', obj.k);
            distance = zeros(1, obj.n);
            for c=1:obj.k
                distance(:,(L==c)') = sum(bsxfun(@minus, obj.matrix((L==c)',:), C(:,c)').^2,2);
            end
            distance = distance';
            s_p = zeros(obj.n, 1);
            for c=1:obj.k
                s_p((L==c),:) = 1/size(L==c,2) + distance((L==c)')./sum(distance);
            end
        end

        function [part, centers] = kmeans(obj, P, k, w)
            if obj.isFastKmeans
                [part, t_centers] = kmeanspp_weighted(P', k, w);
                centers = t_centers';
                %[part, centers, ~] = fkmeans(P, k);
            else
              [part, centers, ~, ~] = Ckmeans(P, k, w,'distance', 'sqeuclidean', ...
                      'maxiter', obj.maxiter, 'emptyaction', 'singleton', 'display', obj.display, 'onlinephase', 'off');

            end
        end

        function [part, centers] = kmeansOneIter(obj, P, k, w)
            if obj.isFastKmeans
                [part, t_centers] = kmeanspp_weighted(P', k, w);
                centers = t_centers';
            else
              [part, centers, ~, ~] = Ckmeans(P, k, w,'distance', 'sqeuclidean', ...
                      'maxiter', 1, 'emptyaction', 'singleton', 'display', obj.display, 'onlinephase', 'off');

            end
        end


        function [error, energy] = getError(obj, subSample, weights)
            [~, ccenter] = obj.kmeans(subSample, obj.k, weights);
            energy = obj.computeEnergy(ccenter);
            error = energy/obj.optEnergy - 1;
        end

        function [error, energy] = nonUniformCoreset(obj)
            try
                weights = 1./(obj.sensitivity()*obj.t);
                sample_idx = randsample(obj.n, obj.t, true, obj.sensitivity());
                coreset = obj.matrix(sample_idx,:);
                [error, energy] = obj.getError(coreset, weights(sample_idx));
            catch err
                disp(err);
                disp(err.stack);
                disp(err.stack.file);
                disp(err.stack.name);
                disp(err.stack.line);
                error = 0;
                energy = 0;
            end
        end

        function [error, energy] = uniformCoreset(obj)
            try
                coreset = datasample(obj.matrix, obj.t);
                [error, energy] = obj.getError(coreset, ones(obj.t, 1));
            catch err
                disp(err);
                disp(err.stack);
                disp(err.stack.file);
                disp(err.stack.name);
                disp(err.stack.line);
                error = 0;
                energy = 0;
            end
        end

        function [error, energy] = kmeanCoreset(obj)
            try
               [m_part, t_centers] = obj.kmeansOneIter(obj.matrix, obj.t, ones(obj.n, 1));

                weights = zeros(obj.t, 1);
                for i=1:obj.t
                    weights(i) = sum(m_part==i);
                end
                [error, energy] = obj.getError(t_centers, weights);
            catch err
                disp(err);
                disp(err.stack);
                err.stack.file
                disp(err.stack.name);
                disp(err.stack.line);
                error = 0;
                energy = 0;
            end
        end

        function obj = run(obj)
            if (obj.synthetic)
                obj.matrix = rand(obj.n, obj.d);
            else
                obj.n = size(obj.matrix, 1);
                obj.d = size(obj.matrix, 2);
            end

            % Compute error based on coreset techniques
            tic;
            [obj.uniformError, obj.uniformEnergy]            = obj.uniformCoreset();
            obj.uniformTime = toc;
            tic;
            [obj.nonUniformError, obj.nonUniformEnergy]      = obj.nonUniformCoreset();
            obj.nonUniformTime = toc;
            tic;
            [obj.kmeansError, obj.kmeansEnergy]              = obj.kmeanCoreset();
            obj.kmeanCoresetTime = toc;
        end
    end

    methods (Static)

        function csn_baja()
            load csn_feature_matrix;
            fileName = '/users/c0rwin/csn_baja.csv';
            for clusters=15:5:25
                test =  KmeansCoresetTest;
                mat = testingMat';
                maxi = 10;

                [~, optCenters] = test.kmeans(mat, clusters, ones(size(mat, 1), 1));
                opt = computeEnergy(mat, optCenters, clusters);

                test.matrix = mat;
                test.fileName = fileName;
                test.toExcel = false;
%                 test.setTestField('caseNo', 1:10);
                test.setTestField('synthetic', false);
                test.setTestField('maxiter', maxi);
                test.setTestField('k', clusters);
                test.setTestField('t', [50:10:100 250:250:5000]);
                test.setTestField('optEnergy', opt);
%                 test.reportFields = { 'caseNo', 'n', 'd', 'k', 't', 'optEnergy', 'uniformEnergy', 'nonUniformEnergy', 'kmeansEnergy', ...
%                     'uniformError', 'nonUniformError', 'kmeansError'};
                test.reportFields = { 'n', 'd', 'k', 't', 'optEnergy', 'uniformEnergy', 'nonUniformEnergy', 'kmeansEnergy', ...
                    'uniformError', 'nonUniformError', 'kmeansError'};

                test.runCartesianProduct( );
            end
        end

        function reuters911()
            load Reuters911;
            for clusters=[3:15 20:5:25]
                test =  KmeansCoresetTest;
                fileName = ['/users/c0rwin/reuters_' num2str(clusters) '.csv'];
                mat = Problem.A';
                maxi = 1;

                [~, optCenters] = test.kmeans(mat, clusters, ones(size(mat, 1), 1));
                opt = computeEnergy(mat, optCenters, clusters);

                test.matrix = mat;
                test.fileName = fileName;
                test.toExcel = false;
                test.setTestField('caseNo', 1:10);
                test.setTestField('synthetic', false);
                test.setTestField('maxiter', maxi);
                test.setTestField('k', clusters);
                test.setTestField('t', 250:250:2500);
                test.setTestField('optEnergy', opt);
                test.reportFields = { 'caseNo', 'n', 'd', 'k', 't', 'optEnergy', 'uniformEnergy', 'nonUniformEnergy', 'kmeansEnergy', ...
                    'uniformError', 'nonUniformError', 'kmeansError'};

                test.runCartesianProduct( );
            end
        end

       function pendigits()
            load data/pendigits.mat;
            fileName = '/users/c0rwin/pendigits_k25_t200.csv';
            for maxi=5:5:100
                clusters = 25;
                test =  KmeansCoresetTest;
                test.maxiter = maxi;
                mat = double([Xte'; Xtr']);

                tic;
                [~, optCenters] = test.kmeans(mat, clusters, ones(size(mat, 1), 1));
                kmeanTime = toc;
                opt = computeEnergy(mat, optCenters, clusters);

                test.matrix = mat;
                test.fileName = fileName;
                test.toExcel = false;

                test.setTestField('caseNo', 1);
                test.setTestField('synthetic', false);
                test.setTestField('maxiter', maxi);
                test.setTestField('k', clusters);
                test.setTestField('t', 200);
                test.setTestField('optEnergy', opt);
                test.setTestField('fullKmeanTime', kmeanTime);
                test.reportFields = { 'caseNo', 'n', 'd', 'k', 't', 'maxiter', 'optEnergy', 'uniformEnergy', 'nonUniformEnergy', 'kmeansEnergy', ...
                    'uniformError', 'nonUniformError', 'kmeansError', 'fullKmeanTime', 'uniformTime', 'nonUniformTime', 'kmeanCoresetTime'};

                test.runCartesianProduct( );
            end
       end

        function mnist()
            load mnist_all;
            fileName = '/users/c0rwin/mnist_final.csv';
            for clusters=15:5:25
                for maxi=10:10:100
                    % clusters = 25;
                    test =  KmeansCoresetTest;
                    test.maxiter = maxi;
                    % test.isFastKmeans = true;
                    mat = double([train1; train2; train3; train4; train5; train6; train7; train8; train9; train0]);

                    tic;
                    [~, optCenters] = test.kmeans(mat, clusters, ones(size(mat, 1), 1));
                    kmeanTime = toc;
                    opt = computeEnergy(mat, optCenters, clusters);

                    test.matrix = mat;
                    test.fileName = fileName;
                    test.toExcel = false;

                    test.setTestField('caseNo', 1);
                    test.setTestField('synthetic', false);
                    test.setTestField('maxiter', maxi);
                    test.setTestField('k', clusters);
                    test.setTestField('t', 500);
                    test.setTestField('optEnergy', opt);
                    test.setTestField('fullKmeanTime', kmeanTime);
                    test.reportFields = { 'caseNo', 'n', 'd', 'k', 't', 'maxiter', 'optEnergy', 'uniformEnergy', 'nonUniformEnergy', 'kmeansEnergy', ...
                        'uniformError', 'nonUniformError', 'kmeansError', 'fullKmeanTime', 'uniformTime', 'nonUniformTime', 'kmeanCoresetTime'};

                    test.runCartesianProduct( );
                end
            end
        end

       function nips()
            load nips_1-17;
            fileName = '/users/c0rwin/nips_k25_t200.csv';
            for maxi=5:5:100
                clusters = 25;
                test =  KmeansCoresetTest;
                mat = aw_counts';
                test.maxiter = maxi;

                tic;
                [~, optCenters] = test.kmeans(mat, clusters, ones(size(mat, 1), 1));
                kmeanTime = toc;
                opt = computeEnergy(mat, optCenters, clusters);

                test.matrix = mat;
                test.fileName = fileName;
                test.toExcel = false;

                test.setTestField('caseNo', 1);
                test.setTestField('synthetic', false);
                test.setTestField('maxiter', maxi);
                test.setTestField('k', clusters);
                test.setTestField('t', 200);
                test.setTestField('optEnergy', opt);
                test.setTestField('fullKmeanTime', kmeanTime);
                test.reportFields = { 'caseNo', 'n', 'd', 'k', 't', 'maxiter', 'optEnergy', 'uniformEnergy', 'nonUniformEnergy', 'kmeansEnergy', ...
                    'uniformError', 'nonUniformError', 'kmeansError', 'fullKmeanTime', 'uniformTime', 'nonUniformTime', 'kmeanCoresetTime'};

                test.runCartesianProduct( );
            end
        end
    end
end
