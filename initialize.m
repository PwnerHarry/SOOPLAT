function initialize()
rng('shuffle');
main_path = fileparts(mfilename('fullpath'));
cd(main_path);
addpath(genpath(fullfile(main_path, 'Algorithms')));
addpath(genpath(fullfile(main_path, 'Public')));
addpath(genpath(fullfile(main_path, 'Benchmarks')));
addpath(genpath(fullfile(main_path, 'Accessories')));
addpath(genpath(fullfile(main_path, 'GUI')));
end