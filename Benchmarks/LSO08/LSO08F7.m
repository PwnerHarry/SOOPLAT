classdef LSO08F7 < PROBLEM
    properties
        ff
    end
    methods
        function obj = LSO08F7(dimension)
            javaclasspath('Benchmarks/LSO08/FractalFunctions.jar');
            % clear java;
            obj.name = 'LSO08F7';
            obj.idealfitness = 0;
            load('Benchmarks/LSO08/LSO08F7.mat');
            obj.idealsolution = NaN(1, dimension);% TODO: make sure they are consistent
            obj.ff = FastFractal('DoubleDip', 3, 1, 1, dimension);
            obj.dimension = dimension;
            obj.lowerbound = -1 * ones(1, dimension);
            obj.upperbound = 1 * ones(1, dimension);
            obj.functionhandle = @(x)obj.ff.evaluate(x);
            obj.idealgroups = {1: dimension};
        end
    end
end