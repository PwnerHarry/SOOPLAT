classdef LSO08F2 < PROBLEM
    methods
        function obj = LSO08F2(dimension)
            obj.name = 'LSO08F2';
            load('Benchmarks/LSO08/LSO08F2.mat');
            if length(o) >= dimension
                obj.idealsolution = o(1: dimension);
            else
                obj.idealsolution = -100 + 200 * rand(1, dimension);
            end
            obj.dimension = dimension;
            obj.lowerbound = -100 * ones(1, dimension);
            obj.upperbound = 100 * ones(1, dimension);
            obj.functionhandle = @(x)obj.schwefel_func(obj.shift(x));
            obj.idealgroups = {1: dimension};
        end
        function fit = schwefel_func(~, x)
            fit = schwefel(x);
        end
    end
end