classdef LSO10F19 < LSO10
    methods
        function obj = LSO10F19(dimension)
            obj.name = 'LSO10F19';
            load('Benchmarks/LSO10/LSO10F19.mat');
            obj.idealsolution = o(1: dimension);
            if dimension ~= 1000
                warning('LSO10F19 supports only D = 1000, automatically corrected');
                dimension = 1000;
            end
            obj.dimension = dimension;
            obj.lowerbound = -100 * ones(1, dimension);
            obj.upperbound = 100 * ones(1, dimension);
            obj.functionhandle = @(x)obj.schwefel_func(obj.shift(x));
            obj.idealgroups = {1: dimension};
            obj.idealfitness = geteps(obj);
        end
    end
end