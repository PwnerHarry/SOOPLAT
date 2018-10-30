function T = getBenchmarkResults(varargin)
% T = getBenchmarkResults('DECCDG2', 'LSO08', 10000, '../SOOPLAT-BenchmarkResults/LSO08/D10000/DECCDG2');
if nargin == 3
    algorithm = varargin{1};
    suite = varargin{2};
    dimension = varargin{3};
    folder = './Data';
elseif nargin == 4
    algorithm = varargin{1};
    suite = varargin{2};
    dimension = varargin{3};
    folder = varargin{4};
end
if strcmp(suite, 'LSO13')
    T = NaN(15, 2);
	if dimension ~= 1000
        warning('LSO13 supports normally D1000, automatically corrected');
    end
    for i = 1: 15
        problem = sprintf('%sF%d', suite, i);
        if i == 13
            dimension = 905;
        else
            dimension = 1000;
        end
        [~, temp] = collectData(51, folder, algorithm, problem, dimension, 3e6);
        results = temp(:, end);
        T(i, 1) = mean(results);
        T(i, 2) = std(results);
    end
elseif strcmp(suite, 'LSO10')
    T = NaN(20, 2);
    if dimension ~= 1000
        warning('LSO10 supports D1000, automatically corrected');
        dimension = 1000;
    end
    for i = 1: 20
        problem = sprintf('%sF%d', suite, i);
        
        [~, temp] = collectData(51, folder, algorithm, problem, dimension, 56);
        results = temp(:, end);
        T(i, 1) = mean(results);
        T(i, 2) = std(results);
    end
elseif  strcmp(suite, 'LSO08')
    T = NaN(7, 2);
    for i = 1: 7
        problem = sprintf('%sF%d', suite, i);
        [~, temp] = collectData(51, folder, algorithm, problem, dimension, 56);
        results = temp(:, end);
        T(i, 1) = mean(results);
        T(i, 2) = std(results);
    end
end
end

