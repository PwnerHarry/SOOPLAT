classdef LSO13F12 < LSO13
    methods
        function obj = LSO13F12(dimension)
            obj.name = 'LSO13F12';
            load('Benchmarks/LSO13/LSO13F12.mat');
            obj.idealsolution = xopt';
            if dimension ~= 1000
                warning('LSO13F12 supports only D = 1000, automatically corrected');
                dimension = 1000;
            end
            obj.dimension = dimension;
            obj.lowerbound = -100 * ones(1, dimension);
            obj.upperbound = 100 * ones(1, dimension);
            obj.functionhandle = @(x)obj.rosenbrock(obj.shift(x)' + 1);
            obj.idealgroups = {1: dimension};
            obj.idealfitness = geteps(obj);
        end
    end
end