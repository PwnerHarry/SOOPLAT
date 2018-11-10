classdef LSO10F13 < LSO10
    properties
        p
    end
    methods
        function obj = LSO10F13(dimension)
            obj.name = 'LSO10F13';
            load('Benchmarks/LSO10/LSO10F13.mat');
            obj.idealsolution = o(1: dimension);
            if dimension ~= 1000
                warning('LSO10F13 supports only D = 1000, automatically corrected');
                dimension = 1000;
            end
            obj.p = p;
            obj.dimension = dimension;
            obj.lowerbound = -100 * ones(1, dimension);
            obj.upperbound = 100 * ones(1, dimension);
            obj.functionhandle = @(x)obj.rosenbrock_group10_func(obj.shift(x));
            obj.idealgroups = {};
            for k = 1: 10
                index = (50 * (k - 1) + 1): 50 * k;
                obj.idealgroups{k} = obj.p(index);
            end
            obj.idealseparables = obj.p(501: end);
            obj.idealfitness = geteps(obj);
        end
        function fit = rosenbrock_group10_func(obj, x)
            fit = 0;
            for k = 1: 10
                index = (50 * (k - 1) + 1): 50 * k;
                fit = fit + obj.rosenbrock_func(x(:, obj.p(index)));
            end
            fit = fit + obj.sphere_func(x(:, obj.p(501: end)));
        end
    end
end