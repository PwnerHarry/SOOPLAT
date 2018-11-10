classdef LSO10F6 < LSO10
    properties
        p
        M
    end
    methods
        function obj = LSO10F6(dimension)
            obj.name = 'LSO10F6';
            load('Benchmarks/LSO10/LSO10F6.mat');
            obj.idealsolution = o(1: dimension);
            if dimension ~= 1000
                warning('LSO10F6 supports only D = 1000, automatically corrected');
                dimension = 1000;
            end
            obj.p = p;
            obj.M = M;
            obj.dimension = dimension;
            obj.lowerbound = -32 * ones(1, dimension);
            obj.upperbound = 32 * ones(1, dimension);
            obj.functionhandle = @(x)obj.ackley_group1_rot_func(obj.shift(x));
            obj.idealgroups = {obj.p(1: 50)};
            obj.idealseparables = obj.p(51: end);
            obj.idealfitness = geteps(obj);
        end
        function fit = ackley_group1_rot_func(obj, x)
            fit = 1e6 * obj.ackley_func(x(:, obj.p(1: 50)) * obj.M) + obj.ackley_func(x(:, obj.p(51: end)));
        end
    end
end