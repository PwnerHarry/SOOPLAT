classdef LSO10F7 < LSO10
    properties
        p
    end
    methods
        function obj = LSO10F7(dimension)
            obj.name = 'LSO10F7';
            load('Benchmarks/LSO10/LSO10F7.mat');
            obj.idealsolution = o(1: dimension);
            if dimension ~= 1000
                warning('LSO10F7 supports only D = 1000, automatically corrected');
                dimension = 1000;
            end
            obj.p = p;
            obj.dimension = dimension;
            obj.lowerbound = -100 * ones(1, dimension);
            obj.upperbound = 100 * ones(1, dimension);
            obj.functionhandle = @(x)obj.schwefel_group1_func(obj.shift(x));
            obj.idealgroups = {obj.p(1: 50)};
            obj.idealseparables = obj.p(51: end);
            obj.idealfitness = geteps(obj);
        end
        function fit = schwefel_group1_func(obj, x)
            fit = 1e6 * obj.schwefel_func(x(:, obj.p(1: 50))) + obj.sphere_func(x(:, obj.p(51: end)));
        end
    end
end
