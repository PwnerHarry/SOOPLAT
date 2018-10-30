classdef LSO10F15 < LSO10
    properties
        p
        M
    end
    methods
        function obj = LSO10F15(dimension)
            obj.name = 'LSO10F15';
            load('Benchmarks/LSO10/LSO10F15.mat');
            obj.idealsolution = o(1: dimension);
            if dimension ~= 1000
                warning('LSO10F15 supports only D = 1000, automatically corrected');
                dimension = 1000;
            end
            obj.p = p;
            obj.M = M;
            obj.dimension = dimension;
            obj.lowerbound = -5 * ones(1, dimension);
            obj.upperbound = 5 * ones(1, dimension);
            obj.functionhandle = @(x)obj.rastrigin_group20_rot_func(obj.shift(x));
            obj.idealgroups = {};
            for k = 1: 20
                index = (50 * (k - 1) + 1): 50 * k;
                obj.idealgroups{k} = obj.p(index);
            end
        end
        function fit = rastrigin_group20_rot_func(obj, x)
            fit = 0;
            for k = 1: 20
                index = (50 * (k - 1) + 1): 50 * k;
                fit = fit + obj.rastrigin_func(x(:, obj.p(index)) * obj.M);
            end
        end
    end
end