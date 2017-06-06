classdef streamingCoresetTest < Test

    properties
        caseNo;

        stream;

        k;

        n;

        d;

        t;

        leafSize;

        kmeansError;

        kmeansEnergy;

        uniformError;

        uniformEnergy;

        nonUniformError;

        nonUniformEnergy;

        optEnergy;

        matrix;

        beta;

        partition;

        isFastKmeans = true;
    end

    methods
        function obj = streamingCoresetTest(k)
            obj.k = k;
        end

        function [part, centers] = kmeans(obj, P, k, w)
            if obj.isFastKmeans
                %[part, t_centers] = kmeanspp_weighted(P', k, w);
                [t_centers, part] = kmeanspp_fast_weighted(P', k, w');
                centers = t_centers';
                %[part, centers, ~] = fkmeans(P, k);
            else
              [part, centers, ~, ~] = Ckmeans(P, k, w,'distance', 'sqeuclidean', ...
                      'maxiter', 100, 'emptyaction', 'singleton', 'display', 'off', 'onlinephase', 'off');

            end
        end

         function energy = computeEnergy(obj, centers)
             dist = zeros(size(obj.matrix, 1), obj.k);
             for c=1:obj.k
                 dist(:,c) = sum(bsxfun(@minus, obj.matrix, centers(c,:)).^2, 2);
             end
             energy = sum(min(dist, [], 2));
         end


        function [error, energy] = getError(obj, subSample, weights)
            try
                rows = size(subSample, 1);
                [~, ccenter] = obj.kmeans(subSample, obj.k, weights(1:rows,:));
                energy = obj.computeEnergy(ccenter);
                error = energy/obj.optEnergy - 1;
            catch err
                disp(err);
                energy = -1;
                error = -1;
            end
        end

        function [error, energy] = kmeanCoreset(obj)
            try
                obj.stream = Stream();
                obj.stream.coresetAlg = kmeansCoreset(obj.t, 10);
                obj.stream.leafSize = obj.leafSize;
                [result] = obj.stream.computeCoreset(PointFunctionSet(obj.matrix));
                weights = result.W.m;
                t_centers = result.M.m;
                [error, energy] = obj.getError(t_centers, weights);
            catch err
                disp(err);
                error = 0;
                energy = 0;
            end
        end

        function [error, energy] = uniformCoreset(obj)
            try
                obj.stream = Stream();
                obj.stream.coresetAlg = uniformCoreset(obj.t);
                obj.stream.leafSize = obj.leafSize;
                [result] = obj.stream.computeCoreset(PointFunctionSet(obj.matrix));
                weights = result.W.m;
                t_centers = result.M.m;
                [error, energy] = obj.getError(t_centers, weights);
            catch err
                disp(err);
                error = 0;
                energy = 0;
            end
        end

        function [error, energy] = nonUniformCoreset(obj)
            try
                obj.stream = Stream();
                obj.stream.coresetAlg = nonUniformCoreset(obj.k, obj.t, 100);
                obj.stream.leafSize = obj.leafSize;
                [result] = obj.stream.computeCoreset(PointFunctionSet(obj.matrix));
                weights = result.W.m;
                t_centers = result.M.m;
                [error, energy] = obj.getError(t_centers, weights);
            catch err
                disp(err);
                error = 0;
                energy = 0;
            end
        end

        function [error, energy] = dannyCoreset(obj)
            try
                obj.stream = Stream();
                dannyCoreset = KMedianCoresetAlg();
                dannyCoreset.coresetType = KMedianCoresetAlg.quadraticInK;
                dannyCoreset.bicriteriaAlg.robustAlg.beta = obj.beta;
                dannyCoreset.bicriteriaAlg.robustAlg.partitionFraction = obj.partition;
                dannyCoreset.bicriteriaAlg.robustAlg.costMethod = ClusterVector.sumSqDistanceCost;
                dannyCoreset.t = obj.t;
                dannyCoreset.k = obj.k;

%                 obj.stream.coresetAlg = nonUniformCoreset(obj.k, obj.t, 100);
                obj.stream.coresetAlg = dannyCoreset;
                obj.stream.leafSize = obj.leafSize;
                [result] = obj.stream.computeCoreset(PointFunctionSet(obj.matrix));
                weights = result.W.m;
                t_centers = result.M.m;
                [error, energy] = obj.getError(t_centers, weights);
            catch err
                disp(err);
                error = 0;
                energy = 0;
            end
        end

        function obj = run(obj)
            obj.n = size(obj.matrix, 1);
            obj.d = size(obj.matrix, 2);

            [obj.nonUniformError, obj.nonUniformEnergy] = obj.dannyCoreset();
            try
                [obj.kmeansError, obj.kmeansEnergy] = obj.kmeanCoreset();
            catch err
                disp(err);
            end
            [obj.uniformError, obj.uniformEnergy] = obj.uniformCoreset();
        end
    end

    methods (Static)

        function random_streaming()
            fileName = '/users/bartem/random_streaming.csv';
            mat = rand(100000, 100);

            clusters = 20;

            test = streamingCoresetTest(clusters);

            [~, optCenters] = test.kmeans(mat, clusters, ones(size(mat, 1), 1));
            opt = computeEnergy(mat, optCenters, clusters);

            test.fileName = fileName;
            test.toExcel = false;
            test.matrix = mat;

            test.setTestField('caseNo', 1:3);
            test.setTestField('leafSize', 1000);
            test.setTestField('t', 50:25:500);
            test.setTestField('k', clusters);
            test.setTestField('beta', 20);
            test.setTestField('partition', 1/2);
            test.setTestField('optEnergy', opt);

            test.reportFields = { 'caseNo', 'n', 'd', 'k', 't', 'leafSize', 'beta', 'partition', 'kmeansError', 'uniformError', 'nonUniformError'};

            test.runCartesianProduct( );
        end

        function mnist_streaming()
            load 'data/mnist_all';
            fileName = '/users/bartem/mnist.csv';
            mat = double([train1; train2; train3; train4; train5; train6; train7; train8; train9; train0]);

            for clusters = 10:5:25;

                test = streamingCoresetTest(clusters);

                %[~, optCenters] = test.kmeans(mat, clusters, ones(size(mat, 1), 1));
                [~, optCenters, ~, ~] = Ckmeans(mat, clusters, ones(size(mat, 1), 1), 'distance', 'sqeuclidean', ...
                      'maxiter', 100, 'emptyaction', 'singleton', 'display', 'off', 'onlinephase', 'off');
                opt = computeEnergy(mat, optCenters, clusters);

                test.fileName = fileName;
                test.toExcel = false;
                test.matrix = mat;

                test.setTestField('caseNo', 1:10);
                test.setTestField('leafSize', [250;500;1000;2000;5000]);
                test.setTestField('t', 50:25:500);
                test.setTestField('k', clusters);
                test.setTestField('beta', [10;15;20;25]);
                test.setTestField('partition', 1/2);
                test.setTestField('optEnergy', opt);

                test.reportFields = { 'caseNo', 'n', 'd', 'k', 't', 'leafSize', 'beta', 'partition', 'kmeansError', 'uniformError', 'nonUniformError'};

                test.runCartesianProduct( );
            end
        end

        function pendigits_streaming()
            fileName = '/users/c0rwin/pendigits_streamingSmall.csv';
            load data/pendigits.mat;

            mat = double([Xte'; Xtr']);

            for clusters = 10:5:25;

                test = streamingCoresetTest(clusters);

                [~, optCenters] = test.kmeans(mat, clusters, ones(size(mat, 1), 1));
                opt = computeEnergy(mat, optCenters, clusters);

                test.fileName = fileName;
                test.toExcel = false;
                test.matrix = mat;

                test.setTestField('caseNo', 1:5);
                test.setTestField('leafSize', 3500);
                test.setTestField('t', 50:25:500);
                test.setTestField('k', clusters);
                test.setTestField('beta', [10 15 20]);
                test.setTestField('partition', 1/2);
                test.setTestField('optEnergy', opt);

                test.reportFields = { 'caseNo', 'n', 'd', 'k', 't', 'leafSize', 'beta', 'partition', 'kmeansError', 'uniformError', 'nonUniformError'};

                test.runCartesianProduct( );
            end
        end

    function erdos_streaming()
            fileName = '/users/bartem/erdos1.csv';
            load data/Erdos982;

            mat = full(Problem.A);

            for clusters = 10:5:25;

                test = streamingCoresetTest(clusters);

                [~, optCenters] = test.kmeans(mat, clusters, ones(size(mat, 1), 1));
                opt = computeEnergy(mat, optCenters, clusters);

                test.fileName = fileName;
                test.toExcel = false;
                test.matrix = mat;
                test.isFastKmeans =true;

                test.setTestField('caseNo', 1:5);
                test.setTestField('leafSize', 3500);
                test.setTestField('t', 50:25:500);
                test.setTestField('k', clusters);
                test.setTestField('beta', [10 15 20]);
                test.setTestField('partition', 1/2);
                test.setTestField('optEnergy', opt);

                test.reportFields = { 'caseNo', 'n', 'd', 'k', 't', 'leafSize', 'beta', 'partition', 'kmeansError', 'uniformError', 'nonUniformError'};

                test.runCartesianProduct( );
            end
    end

            function reuters_streaming()
            fileName = '/users/bartem/reuters.csv';
            load data/Reuters911;

            mat = full(Problem.A);

            for clusters = 10:5:25;

                test = streamingCoresetTest(clusters);

                [~, optCenters] = test.kmeans(mat, clusters, ones(size(mat, 1), 1));
                opt = computeEnergy(mat, optCenters, clusters);

                test.fileName = fileName;
                test.toExcel = false;
                test.matrix = mat;
                test.isFastKmeans =true;

                test.setTestField('caseNo', 1:5);
                test.setTestField('leafSize', 3500);
                test.setTestField('t', 50:25:500);
                test.setTestField('k', clusters);
                test.setTestField('beta', 20);
                test.setTestField('partition', 1/2);
                test.setTestField('optEnergy', opt);

                test.reportFields = { 'caseNo', 'n', 'd', 'k', 't', 'leafSize', 'beta', 'partition', 'kmeansError', 'uniformError', 'nonUniformError'};

                test.runCartesianProduct( );
            end
            end

        function imdb_streaming()
            fileName = '/users/bartem/imdb.csv';
            load data/web-Stanford;

            mat = Problem.A;

            for clusters = 10:5:25;

                test = streamingCoresetTest(clusters);

                [~, optCenters] = test.kmeans(mat, clusters, ones(size(mat, 1), 1));
                opt = computeEnergy(mat, optCenters, clusters);

                test.fileName = fileName;
                test.toExcel = false;
                test.matrix = mat;
                test.isFastKmeans =true;

                test.setTestField('caseNo', 1:5);
                test.setTestField('leafSize', 3500);
                test.setTestField('t', 50:25:500);
                test.setTestField('k', clusters);
                test.setTestField('beta', 20);
                test.setTestField('partition', 1/2);
                test.setTestField('optEnergy', opt);

                test.reportFields = { 'caseNo', 'n', 'd', 'k', 't', 'leafSize', 'beta', 'partition', 'kmeansError', 'uniformError', 'nonUniformError'};

                test.runCartesianProduct( );
            end
        end


        function imdb1_streaming()
            fileName = '/users/bartem/imdb.csv';
            load 'data/IMDB';

            mat = Problem.A;

            for clusters = 10:5:25;

                test = streamingCoresetTest(clusters);

                [~, optCenters] = test.kmeans(mat, clusters, ones(size(mat, 1), 1));
                opt = computeEnergy(mat, optCenters, clusters);

                test.fileName = fileName;
                test.toExcel = false;
                test.matrix = mat;

                test.setTestField('caseNo', 1);
                test.setTestField('leafSize', 500);
                test.setTestField('t', 50:25:500);
                test.setTestField('k', clusters);
                test.setTestField('beta', [10 15 20]);
                test.setTestField('partition', 1/2);
                test.setTestField('optEnergy', opt);

                test.reportFields = { 'caseNo', 'n', 'd', 'k', 't', 'leafSize', 'beta', 'partition', 'kmeansError', 'uniformError', 'nonUniformError'};

                test.runCartesianProduct( );
            end
        end

    end
end
