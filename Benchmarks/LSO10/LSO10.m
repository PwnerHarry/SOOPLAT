classdef LSO10 < PROBLEM
    methods
        function obj = LSO10()
        end
        function fit = elliptic_func(~, x)
            D = size(x, 2);
            fit = 0;
            for i = 1: D
                fit = fit + 1e+6 .^ ((i - 1) / (D - 1)) .* x(:, i) .^ 2;
            end
        end
        function fit = rastrigin_func(~, x)
            fit = sum(x .* x - 10 * cos(2 * pi * x) + 10, 2);
        end
        function fit = ackley_func(~, x)
            D = size(x, 2);
            fit = sum(x .* x, 2);
            fit = 20 * (1 - exp(-0.2 .* sqrt(fit ./ D))) + exp(1) - exp(sum(cos(2 .* pi .* x), 2) ./ D);
        end
        function fit = schwefel_func(~, x)
            [ps, D] = size(x);
            fit = 0;
            S = zeros(ps, 1);
            for i = 1: D
                S = S + x(:, i);
                fit = fit + S .^ 2;
            end
        end
        function fit = rosenbrock_func(~, x)
            D = size(x, 2);
            fit = sum(100 .* (x(:, 1: D - 1) .^ 2 - x(:, 2: D)) .^ 2 + (x(:, 1: D - 1) - 1) .^ 2, 2);
        end
        function fit = sphere_func(~, x)
            fit = sum(x .* x, 2);
        end
    end
end