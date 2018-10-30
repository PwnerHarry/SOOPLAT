classdef LSO13F3 < LSO13
    methods
        function obj = LSO13F3(dimension)
            obj.name = 'LSO13F3';
            load('Benchmarks/LSO13/LSO13F3.mat');
            obj.idealsolution = xopt';
            if dimension ~= 1000
                warning('LSO13F3 supports only D = 1000, automatically corrected');
                dimension = 1000;
            end
            obj.dimension = dimension;
            obj.lowerbound = -32 * ones(1, dimension);
            obj.upperbound = 32 * ones(1, dimension);
            obj.functionhandle = @(x)obj.ackley(obj.shift(x)');
            obj.idealseparables = 1: dimension;
        end
    end
end