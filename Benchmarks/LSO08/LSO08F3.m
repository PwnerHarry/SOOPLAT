classdef LSO08F3 < PROBLEM
    methods
        function obj = LSO08F3(dimension)
            obj.name = 'LSO08F3';
            load('Benchmarks/LSO08/LSO08F3.mat');
            if length(o) >= dimension
                obj.idealsolution = o(1: dimension);
            else
                obj.idealsolution = -90 + 180 * rand(1, dimension);
            end
            obj.dimension = dimension;
            obj.lowerbound = -100 * ones(1, dimension);
            obj.upperbound = 100 * ones(1, dimension);
            obj.functionhandle = @(x)obj.rosenbrock_func(obj.shift(x));
            obj.idealgroups = {};
            obj.idealseparables = 1: dimension;
        end
        function f = rosenbrock_func(~, x)
            % D = size(x, 2);
            % f = sum(100 .* (x(:, 1: D - 1) .^ 2 - x(:, 2: D)) .^ 2 + (x(:, 1: D - 1) - 1) .^ 2, 2);
            f = rosenbrock(x);
        end
    end
end