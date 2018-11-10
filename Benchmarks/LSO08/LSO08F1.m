classdef LSO08F1 < PROBLEM
    methods
        function obj = LSO08F1(dimension)
            obj.name = 'LSO08F1';
            load('Benchmarks/LSO08/LSO08F1.mat');
            if length(o) >= dimension
                obj.idealsolution = o(1: dimension);
            else
                obj.idealsolution = -100 + 200 * rand(1, dimension);
            end
            obj.dimension = dimension;
            obj.lowerbound = -100 * ones(1, dimension);
            obj.upperbound = 100 * ones(1, dimension);
            obj.functionhandle = @(x)obj.sphere_func(obj.shift(x));
            obj.idealgroups = {};
            obj.idealseparables = 1: dimension;
            obj.idealfitness = geteps(obj);
        end
        function fit = sphere_func(~, x)
            fit = sphere(x);
        end
    end
end
