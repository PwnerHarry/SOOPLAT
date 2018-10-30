classdef LSO13 < PROBLEM
    methods
        function obj = LSO13()
        end
        function fit = elliptic(obj, x)
            coefficients = 1e6 .^ linspace(0, 1, size(x, 1));
            X = obj.T_irreg(x);
            fit = coefficients * (X .* X);
        end
        function fit = rastrigin(obj, x)
            x = obj.T_diag(obj.T_asy(obj.T_irreg(x), 0.2), 10);
            fit = 10 * (size(x, 1) - sum(cos(2 * pi * x), 1)) + sum(x .* x, 1);
        end
        function fit = ackley(~, x)
            D = size(x, 1);
            fit = sum(x .* x,1);
            fit = 20 .* (1 - exp(- 0.2 .* sqrt(fit ./ D))) - exp(sum(cos(2 .* pi .* x), 1) ./ D) + exp(1);
        end
        function fit = sphere(~, x)
            fit = sum(x .* x);
        end
        
        function fit = schwefel(obj, x)
            x = obj.T_asy(obj.T_irreg(x), 0.2);
            fit = 0;
            S = zeros(1, size(x, 2));
            for i = 1: size(x, 1)
                S = S + x(i, :);
                fit = fit + S .^ 2;
            end
        end
        function fit = rosenbrock(~, x)
            D = size(x, 1);
            fit = sum(100 .* (x(1: D - 1, :) .^ 2 - x(2: D, :)) .^ 2 + (x(1: D - 1, :) - 1) .^ 2);
        end
        function g = T_asy(~, f, beta)
            [D, popsize] = size(f);
            g = f;
            temp = repmat(beta * linspace(0, 1, D)', 1, popsize);
            ind = f > 0;
            g(ind) = f(ind).^ (1 + temp(ind) .* sqrt(f(ind)));
        end
        
        function g = T_diag(~, f, alpha)
            [D, popsize] = size(f);
            scales = repmat(sqrt(alpha) .^ linspace(0, 1, D)', 1, popsize);
            g = scales .* f;
        end
        
        function g = T_irreg(~, f)
            g = f;
            idx = (f > 0);
            g(idx) = 10 .* log(f(idx));
            g(idx) = exp(0.1 .* (g(idx) + 0.49 * (sin(g(idx)) + sin(0.79 * g(idx)))));
            idx = (f < 0);
            g(idx) = 10 .* log(-f(idx));
            g(idx) = - exp(0.1 .* (g(idx) + 0.49 * (sin(0.55 * g(idx)) + sin(0.31 * g(idx)))));
        end
    end
end