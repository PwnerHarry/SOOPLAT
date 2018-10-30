classdef LSO08F4 < PROBLEM
    methods
        function obj = LSO08F4(dimension)
            obj.name = 'LSO08F4';
            obj.dimension = dimension;
            load('Benchmarks/LSO08/LSO08F4.mat');
            if length(o) >= dimension
                obj.idealsolution = o(1: dimension);
            else
                obj.idealsolution = -5 + 10 * rand(1, dimension);
            end
            obj.lowerbound = -5 * ones(1, dimension);
            obj.upperbound = 5 * ones(1, dimension);
            obj.functionhandle = @(x)obj.rastrigin_func(obj.shift(x));
            obj.idealgroups = {};
            obj.idealseparables = 1: dimension;
        end
        function f = rastrigin_func(~, x)
            % f = sum(x .^ 2 - 10 .* cos(2 .* pi .* x) + 10, 2);
            f = rastrigin(x);
        end
    end
end