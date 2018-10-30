classdef LSO10F18 < LSO10
    properties
        p
    end
    methods
        function obj = LSO10F18(dimension)
            obj.name = 'LSO10F18';
            load('Benchmarks/LSO10/LSO10F18.mat');
            obj.idealsolution = o(1: dimension);
            obj.p = p;
            if dimension ~= 1000
                warning('LSO10F18 supports only D = 1000, automatically corrected');
                dimension = 1000;
            end
            obj.dimension = dimension;
            obj.lowerbound = -100 * ones(1, dimension);
            obj.upperbound = 100 * ones(1, dimension);
            obj.functionhandle = @(x)obj.rosenbrock_group20_func(obj.shift(x));
            obj.idealgroups = {};
            for k = 1: 20
                index = (50 * (k - 1) + 1): 50 * k;
                obj.idealgroups{k} = obj.p(index);
            end
        end
        function fit = rosenbrock_group20_func(obj, x)
            fit = 0;
            for k = 1: 20
                index = (50 * (k - 1) + 1): 50 * k;
                fit = fit + obj.rosenbrock_func(x(:, obj.p(index)));
            end
        end
    end
end