classdef LSO13F13 < LSO13
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
        function obj = LSO13F13(dimension)
            obj.name = 'LSO13F13';
            load('Benchmarks/LSO13/LSO13F13.mat');
            obj.idealsolution = xopt';
            if dimension ~= 905
                warning('LSO13F13 supports only D = 905, automatically corrected');
                dimension = 905;
            end
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
            obj.functionhandle = @(x)obj.f13(obj.shift(x)');
            obj.idealgroups = {};
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
        function fit = f13(obj, x)
            fit = 0;
            for i = 1: length(obj.s)
                if i == 1
                    ldim = 1;
                else
                    ldim = obj.c(i - 1) - ((i - 1) * obj.m) + 1;
                end
                udim = obj.c(i) - ((i - 1) * obj.m);
                if obj.s(i) == 25
                    f = obj.schwefel(obj.R25 * x(obj.p(ldim: udim), :));
                elseif obj.s(i) == 50
                    f = obj.schwefel(obj.R50 * x(obj.p(ldim: udim), :));
                elseif obj.s(i) == 100
                    f = obj.schwefel(obj.R100 * x(obj.p(ldim: udim), :));
                end
                fit = fit + obj.w(i) * f;
            end
        end
    end
end