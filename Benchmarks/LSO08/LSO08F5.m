classdef LSO08F5 < PROBLEM
    methods
        function obj = LSO08F5(dimension)
            obj.name = 'LSO08F5';
            load('Benchmarks/LSO08/LSO08F5.mat');
            if length(o) >= dimension
                obj.idealsolution = o(1: dimension);
            else
                obj.idealsolution = -600 + 1200 * rand(1, dimension);
            end
            obj.dimension = dimension;
            obj.lowerbound = -600 * ones(1, dimension);
            obj.upperbound = 600 * ones(1, dimension);
            obj.functionhandle = @(x)obj.griewank_func(obj.shift(x));
            obj.idealgroups = {1: dimension};
        end
        function f = griewank_func(~, x)
%             D = size(x, 2);
%             f = 1;
%             for i = 1: D
%                 f = f .* cos(x(:, i) ./ sqrt(i));
%             end
%             f = sum(x .* x, 2) ./ 4000 - f + 1;
            f = griewank(x);
        end
    end
end