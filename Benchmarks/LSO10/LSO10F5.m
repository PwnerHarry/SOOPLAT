classdef LSO10F5 < LSO10
    properties
        p
        M
    end
    methods
        function obj = LSO10F5(dimension)
            obj.name = 'LSO10F5';
            load('Benchmarks/LSO10/LSO10F5.mat');
            obj.idealsolution = o(1: dimension);
            if dimension ~= 1000
                warning('LSO10F5 supports only D = 1000, automatically corrected');
                dimension = 1000;
            end
            obj.p = p;
            obj.M = M;
            obj.dimension = dimension;
            obj.lowerbound = -5 * ones(1, dimension);
            obj.upperbound = 5 * ones(1, dimension);
            obj.functionhandle = @(x)obj.rastrigin_group1_rot_func(obj.shift(x));
            obj.idealgroups = {obj.p(1: 50)};
            obj.idealseparables = obj.p(51: end);
        end
        function fit = rastrigin_group1_rot_func(obj, x)
            fit = 1e6 * obj.rastrigin_func(x(:,obj.p(1: 50)) * obj.M) + obj.rastrigin_func(x(:, obj.p(51: end)));
        end
    end
end