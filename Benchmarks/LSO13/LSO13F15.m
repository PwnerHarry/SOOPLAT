classdef LSO13F15 < LSO13
    methods
        function obj = LSO13F15(dimension)
            obj.name = 'LSO13F15';
            load('Benchmarks/LSO13/LSO13F15.mat');
            obj.idealsolution = xopt';
            if dimension ~= 1000
                warning('LSO13F15 supports only D = 1000, automatically corrected');
                dimension = 1000;
            end
            obj.dimension = dimension;
            obj.lowerbound = -100 * ones(1, dimension);
            obj.upperbound = 100 * ones(1, dimension);
            obj.functionhandle = @(x)obj.schwefel(obj.shift(x)');
            obj.idealgroups = {1: dimension};
        end
    end
end