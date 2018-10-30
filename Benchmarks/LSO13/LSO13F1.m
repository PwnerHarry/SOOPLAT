classdef LSO13F1 < LSO13
    methods
        function obj = LSO13F1(dimension)
            obj.name = 'LSO13F1';
            load('Benchmarks/LSO13/LSO13F1.mat');
            obj.idealsolution = xopt';
            if dimension ~= 1000
                warning('LSO13F1 supports only D = 1000, automatically corrected');
                dimension = 1000;
            end
            obj.dimension = dimension;
            obj.lowerbound = -100 * ones(1, dimension);
            obj.upperbound = 100 * ones(1, dimension);
            obj.functionhandle = @(x)obj.elliptic(obj.shift(x)');
            obj.idealseparables = 1: dimension;
        end
    end
end