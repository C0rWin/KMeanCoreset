classdef KmeansCoresetTest < Test

    properties

        caseNo; n; d; k; t; realSize; matrix;

        optEnergy; uniformEnergy; nonUniformEnergy; kmeansEnergy;

        uniformError; nonUniformError; kmeansError;

        kmeansErrorWithSearch;

        maxiter = 100; display = 'off';

        synthetic = true;

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

        function [error, energy] = getError(obj, subSample, weights)
%             [~, ccenter] = Ckmeans(subSample, obj.k, weights, 'distance', 'sqeuclidean', ...
%                 'maxiter', obj.maxiter, 'emptyaction', 'singleton', 'display', obj.display, 'onlinephase', 'off');
            [~, ccenter] = kmeanspp_weighted(subSample', obj.k, weights);
            energy = obj.computeEnergy(ccenter');
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
                error = 0;
                energy = 0;
            end
        end

        function [error, energy] = kmeanCoreset(obj)
            try
                 [m_part, t_centers, ~, ~] = Ckmeans(obj.matrix, obj.t, ones(obj.n, 1),'distance', 'sqeuclidean', ...
                     'maxiter', obj.maxiter, 'emptyaction', 'singleton', 'display', obj.display, 'onlinephase', 'off');
%                 [m_part, t_centers] = kmeanspp_weighted(obj.matrix', obj.t, ones(obj.n, 1));

                weights = zeros(obj.t, 1);
                for i=1:obj.t
                    weights(i) = sum(m_part==i);
                end
                [error, energy] = obj.getError(t_centers, weights);
            catch err
                disp(err);
                error = 0;
                energy = 0;
            end
        end

        function [error, energy] = kmeanCoresetWithSearch(obj)
            try
                perf_t = zeros(obj.t - obj.k, 1);
                for i=obj.k:obj.t
                    try
                        [~, ~, dists_i, ~] = Ckmeans(obj.matrix, i, ones(obj.n, 1), 'distance', 'sqeuclidean', ...
                            'maxiter', obj.maxiter, 'emptyaction', 'singleton', 'display', obj.display, 'onlinephase', 'off');
                        [~, ~, dists_i_k_1, ~] = Ckmeans(obj.matrix, i+obj.k-1, ones(obj.n, 1), 'distance', 'sqeuclidean', ...
                            'maxiter', obj.maxiter, 'emptyaction', 'singleton', 'display', obj.display, 'onlinephase', 'off');
                        perf_t = sum(dists_i) - sum(dists_i_k_1);
                    catch err
                        disp(err);
                        continue;
                    end
                end
                [~, obj.realSize] = min(perf_t);
                [m_part, t_centers, ~, ~] = Ckmeans(obj.matrix, obj.realSize, ones(obj.n, 1), 'distance', 'sqeuclidean', 'maxiter', obj.maxiter, ...
                    'emptyaction', 'singleton', 'display', obj.display, 'onlinephase', 'off');
                weights = zeros(obj.realSize, 1);
                for i=1:obj.realSize;
                    weights(i) = sum(m_part==i);
                end
                [error, energy] = obj.getError(t_centers, weights);
            catch err
                disp(err);
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

            % Compute error based on coreset teechniques
            [obj.uniformError, obj.uniformEnergy]            = obj.uniformCoreset();
            [obj.nonUniformError, obj.nonUniformEnergy]      = obj.nonUniformCoreset();
            [obj.kmeansError, obj.kmeansEnergy]              = obj.kmeanCoreset();
%             obj.kmeansErrorWithSearch   = obj.kmeanCoresetWithSearch();
        end
    end

    methods (Static)
        function mnist()
            load datasets/mnist_all;
            for clusters=5:5:25
                test =  KmeansCoresetTest;
                fileName = ['/users/c0rwin/mnist_' num2str(clusters) '.csv'];
                mat = double([train1; train2; train3; train4; train5; train6; train7; train8; train9; train0]);
                maxi = 10;

%                 [~, optCenters] = Ckmeans(mat, clusters, ones(size(mat, 1), 1), 'distance', 'sqeuclidean', ...
%                     'maxiter', maxi, 'emptyaction', 'singleton', 'display', 'off', 'onlinephase', 'off');
                [~, optCenters] = kmeanspp_weighted(mat', clusters, ones(size(mat, 1), 1));
                opt = computeEnergy(mat, optCenters', clusters);

                test.matrix = mat;
                test.fileName = fileName;
                test.toExcel = false;

                test.setTestField('caseNo', 1:10);
                test.setTestField('synthetic', false);
                test.setTestField('maxiter', maxi);
                test.setTestField('k', clusters);
                test.setTestField('t', 50:10:100);
                test.setTestField('optEnergy', opt);

                test.reportFields = { 'caseNo', 'n', 'd', 'k', 't', ...
                    'optEnergy', 'uniformEnergy', 'nonUniformEnergy', 'kmeansEnergy', ...
                    'uniformError', 'nonUniformError', 'kmeansError'};
                test.runCartesianProduct( );
            end
        end

        function csn_baja()
            load csn_baja_feature_matrix;
            for clusters=[3:15 20:5:25]
                test =  KmeansCoresetTest;
                fileName = ['/users/c0rwin/csn_baja_' num2str(clusters) '.csv'];
                mat = testingMat';
                maxi = 10;

%                 [~, optCenters] = Ckmeans(mat, clusters, ones(size(mat, 1), 1), 'distance', 'sqeuclidean', ...
%                     'maxiter', maxi, 'emptyaction', 'singleton', 'display', 'off', 'onlinephase', 'off');
                [~, optCenters] = kmeanspp_weighted(mat', clusters, ones(size(mat, 1), 1));

                opt = computeEnergy(mat, optCenters', clusters);

                test.matrix = mat;
                test.fileName = fileName;
                test.toExcel = false;
                test.setTestField('caseNo', 1:10);
                test.setTestField('synthetic', false);
                test.setTestField('maxiter', maxi);
                test.setTestField('k', clusters);
                test.setTestField('t', 50:5:100);
                test.setTestField('optEnergy', opt);

                test.reportFields = { 'caseNo', 'n', 'd', 'k', 't', ...
                    'optEnergy', 'uniformEnergy', 'nonUniformEnergy', 'kmeansEnergy', ...
                    'uniformError', 'nonUniformError', 'kmeansError'};
                test.runCartesianProduct( );
            end
        end

        function csn()
            load csn_feature_matrix;
%             test.matrix = double([train1; train2; train3; train4; train5; train6; train7; train8; train9; train0]);
            for clusters=[3:15 20:5:30]
                test =  KmeansCoresetTest;
                fileName = ['/users/c0rwin/csn_' num2str(clusters) '.csv'];
                mat = testingMat';
                maxi = 10;

                [~, optCenters] = Ckmeans(mat, clusters, ones(size(mat, 1), 1), 'distance', 'sqeuclidean', ...
                    'maxiter', maxi, 'emptyaction', 'singleton', 'display', 'off', 'onlinephase', 'off');
                opt = computeEnergy(mat, optCenters, clusters);

                test.matrix = mat;
                test.fileName = fileName;
                test.toExcel = false;

                test.setTestField('synthetic', false);
                test.setTestField('maxiter', maxi);
                test.setTestField('k', clusters);
                test.setTestField('t', 250:250:5000);
                test.setTestField('optEnergy', opt);

                test.reportFields = { 'n', 'd', 'k', 't', ...
                    'optEnergy', 'uniformEnergy', 'nonUniformEnergy', 'kmeansEnergy', ...
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

                 [~, optCenters] = Ckmeans(mat, clusters, ones(size(mat, 1), 1), 'distance', 'sqeuclidean', ...
                     'maxiter', maxi, 'emptyaction', 'singleton', 'display', 'off', 'onlinephase', 'off');
%                 [~, optCenters] = kmeanspp_weighted(mat', clusters, ones(size(mat, 1), 1));

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

                test.reportFields = { 'caseNo', 'n', 'd', 'k', 't', ...
                    'optEnergy', 'uniformEnergy', 'nonUniformEnergy', 'kmeansEnergy', ...
                    'uniformError', 'nonUniformError', 'kmeansError'};
                test.runCartesianProduct( );
            end
        end


        function run_csns()
            KmeansCoresetTest.csn();
            KmeansCoresetTest.csn_baja();
        end
    end
end
