classdef LSO13F14 < LSO13
    properties
        p
        s
        R25
        R50
        R100
        m
        c
        w
    end
    methods
        function obj = LSO13F14(dimension)
            obj.name = 'LSO13F14';
            global LSO13F14_INITIALIZED;
            LSO13F14_INITIALIZED = false;
            if dimension ~= 1000
                warning('LSO13F14 supports only D = 1000, automatically corrected');
                dimension = 1000;
            end
            load('Benchmarks/LSO13/LSO13F14.mat');
            obj.idealsolution = xopt(1: dimension)';
            %obj.idealsolution = zeros(1, dimension);
            obj.p = p;
            obj.s = s;
            obj.R25 = R25;
            obj.R50 = R50;
            obj.R100 = R100;
            obj.m = m;
            obj.c = cumsum(obj.s);
            obj.w = w;
            obj.dimension = dimension;
            obj.lowerbound = -100 * ones(1, dimension);
            obj.upperbound = 100 * ones(1, dimension);
            obj.functionhandle = @(x)obj.f14(x');
            obj.idealgroups = {};
            obj.idealfitness = eps; % obj.idealfitness = geteps(obj);
            for i = 1: length(obj.s)
                if i == 1
                    ldim = 1;
                else
                    ldim = obj.c(i - 1) - ((i - 1) * obj.m) + 1;
                end
                udim = obj.c(i) - ((i - 1) * obj.m);
                obj.idealgroups{i} = obj.p(ldim: udim);
            end
        end
        function fit = f14(obj, x)
            fit = 0;
            for i = 1: length(obj.s)
                if i == 1
                    ldim = 1;
                    ldimshift = 1;
                else
                    ldim = obj.c(i - 1) - ((i - 1) * obj.m) + 1;
                    ldimshift = obj.c(i - 1) + 1;
                end
                udim = obj.c(i) - ((i - 1) * obj.m);
                udimshift = obj.c(i);
                z = x(obj.p(ldim: udim), :) - repmat(obj.idealsolution(ldimshift: udimshift)', 1, size(x, 2));
                if obj.s(i) == 25
                    f = obj.schwefel(obj.R25 * z);
                elseif obj.s(i) == 50
                    f = obj.schwefel(obj.R50 * z);
                elseif obj.s(i) == 100
                    f = obj.schwefel(obj.R100 * z);
                end
                fit = fit + obj.w(i) * f;
            end
        end
    end
end