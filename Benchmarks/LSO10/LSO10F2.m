classdef LSO10F2 < LSO10
    methods
        function obj = LSO10F2(dimension)
            obj.name = 'LSO10F2';
            load('Benchmarks/LSO10/LSO10F2.mat');
            obj.idealsolution = o(1: dimension);
            obj.dimension = dimension;
            obj.lowerbound = -5 * ones(1, dimension);
            obj.upperbound = 5 * ones(1, dimension);
            obj.functionhandle = @(x)obj.rastrigin_func(obj.shift(x));
            obj.idealseparables = 1: dimension;
        end
        
    end
end


