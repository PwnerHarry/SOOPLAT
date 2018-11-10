classdef LSO10F1 < LSO10
    methods
        function obj = LSO10F1(dimension)
            obj.name = 'LSO10F1';
            load('Benchmarks/LSO10/LSO10F1.mat');
            obj.idealsolution = o(1: dimension);
            obj.dimension = dimension;
            obj.lowerbound = -100 * ones(1, dimension);
            obj.upperbound = 100 * ones(1, dimension);
            obj.functionhandle = @(x)obj.elliptic_func(obj.shift(x));
            obj.idealseparables = 1: dimension;
            obj.idealfitness = geteps(obj);
        end
    end
end