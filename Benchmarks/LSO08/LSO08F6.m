classdef LSO08F6 < PROBLEM
    methods
        function obj = LSO08F6(dimension)
            obj.name = 'LSO08F6';
            load('Benchmarks/LSO08/LSO08F6.mat');
            if length(o) >= dimension
                obj.idealsolution = o(1: dimension);
            else
                obj.idealsolution = -30 + 60 * rand(1, dimension);
            end
            obj.dimension = dimension;
            obj.lowerbound = -32 * ones(1, dimension);
            obj.upperbound = 32 * ones(1, dimension);
            obj.functionhandle = @(x)obj.ackley_func(obj.shift(x));
            obj.idealgroups = {};
            obj.idealseparables = 1: dimension;
            obj.idealfitness = geteps(obj);
        end
        function f = ackley_func(~, x)
            % D = size(x, 2);
            % f = 20 * (1 - exp(-0.2 .* sqrt(sum(x .* x, 2) ./ D))) + exp(1) * (1 - exp(sum(cos(2 .* pi .* x), 2) ./ D - 1));
            f = ackley(x);
        end
    end
end