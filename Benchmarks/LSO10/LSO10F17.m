classdef LSO10F17 < LSO10
    properties
        p
    end
    methods
        function obj = LSO10F17(dimension)
            obj.name = 'LSO10F17';
            load('Benchmarks/LSO10/LSO10F17.mat');
            obj.idealsolution = o(1: dimension);
            if dimension ~= 1000
                warning('LSO10F17 supports only D = 1000, automatically corrected');
                dimension = 1000;
            end
            obj.p = p;
            obj.dimension = dimension;
            obj.lowerbound = -100 * ones(1, dimension);
            obj.upperbound = 100 * ones(1, dimension);
            obj.functionhandle = @(x)obj.schwefel_group20_func(obj.shift(x));
            obj.idealgroups = {};
            for k = 1: 20
                index = (50 * (k - 1) + 1): 50 * k;
                obj.idealgroups{k} = obj.p(index);
            end
        end
        function fit = schwefel_group20_func(obj, x)
            fit = 0;
            for k = 1: 20
                index = (50 * (k - 1) + 1): 50 * k;
                fit = fit + obj.schwefel_func(x(:, obj.p(index)));
            end
        end
    end
end