classdef LSO10F3 < LSO10
    methods
        function obj = LSO10F3(dimension)
            obj.name = 'LSO10F3';
            load('Benchmarks/LSO10/LSO10F3.mat');
            obj.idealsolution = o(1: dimension);
            obj.dimension = dimension;
            obj.lowerbound = -32 * ones(1, dimension);
            obj.upperbound = 32 * ones(1, dimension);
            obj.functionhandle = @(x)obj.ackley_func(obj.shift(x));
            obj.idealseparables = 1: dimension;
        end
    end
end