classdef LSO13F2 < LSO13
    methods
        function obj = LSO13F2(dimension)
            obj.name = 'LSO13F2';
            load('Benchmarks/LSO13/LSO13F2.mat');
            obj.idealsolution = xopt';
            if dimension ~= 1000
                warning('LSO13F2 supports only D = 1000, automatically corrected');
                dimension = 1000;
            end
            obj.dimension = dimension;
            obj.lowerbound = -5 * ones(1, dimension);
            obj.upperbound = 5 * ones(1, dimension);
            obj.functionhandle = @(x)obj.rastrigin(obj.shift(x)');
            obj.idealseparables = 1: dimension;
            obj.idealfitness = geteps(obj);
        end
    end
end