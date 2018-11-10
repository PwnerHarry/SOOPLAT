classdef LSO10F4 < LSO10
    properties
        p
        M
    end
    methods
        function obj = LSO10F4(dimension)
            obj.name = 'LSO10F4';
            load('Benchmarks/LSO10/LSO10F4.mat');
            obj.idealsolution = o(1: dimension);
            obj.p = p;
            obj.M = M;
            if dimension ~= 1000
                warning('LSO10F4 supports only D = 1000, automatically corrected');
                dimension = 1000;
            end
            obj.dimension = dimension;
            obj.lowerbound = -100 * ones(1, dimension);
            obj.upperbound = 100 * ones(1, dimension);
            obj.functionhandle = @(x)obj.elliptic_group1_rot_func(obj.shift(x));
            obj.idealgroups = {obj.p(1: 50)};
            obj.idealseparables = obj.p(51: end);
            obj.idealfitness = geteps(obj);
        end
        function fit = elliptic_group1_rot_func(obj, x)
            fit = 1e6 * obj.elliptic_func(x(:, obj.p(1: 50)) * obj.M) + obj.elliptic_func(x(:, obj.p(51: end)));
        end
    end
end