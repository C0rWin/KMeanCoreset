classdef kmeansStreamingTest < Test

    properties
        allC;

        caseNo;

        no;

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

        memory = 0;
    end

    methods
        function obj = kmeansStreamingTest(k)
            obj.k = k;
            obj.allC = [];
        end

        function [part, centers] = kmeans(obj, P, k, w)
%            if obj.isFastKmeans
%                [part, t_centers] = kmeanspp_weighted(P', k, w);
%                centers = t_centers';
                %[part, centers, ~] = fkmeans(P, k);
%            else
              [part, centers, ~, ~] = Ckmeans(P, k, w,'distance', 'sqeuclidean', ...
                      'maxiter', 100, 'emptyaction', 'singleton', 'display', 'off', 'onlinephase', 'off');

%            end
        end

        function kmeanCoreset(obj)
            obj.no = obj.no + 1;
            obj.stream.addPointSet(obj.matrix);
            obj.memory = size(obj.stream.stack.stack, 2);
        end

        function obj = run(obj)
            obj.n = obj.stream.numPointsStreamed + size(obj.matrix.M.m, 1);
            obj.d = size(obj.matrix.M.m, 2);

            obj.kmeanCoreset();
        end
    end

    methods (Static)

        function mnist_streaming()

            fileName = '/users/c0rwin/memory.csv';
            clusters = 20;

            test = kmeansStreamingTest(clusters);

            test.fileName = fileName;
            test.toExcel = false;
            test.stream = Stream();
            test.stream.leafSize = 500;
            test.stream.coresetAlg = kmeansCoreset(200, 100);

            for i=0:10000

                test.setTestField('t', 200);
               test.setTestField('k', clusters);
                test.setTestField('matrix', PointFunctionSet(rand(1000, 100)));

                test.reportFields = { 'caseNo', 'n', 'd', 'k', 'memory'};

                test.runCartesianProduct( );
            end
        end
    end
end
