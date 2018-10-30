function initSOOPLAT()
    % Please do not use "clear functions" anywhere else
    % for the SOOPLAT relies on the persistent variables
    rng('shuffle');
    main_path = fileparts(mfilename('fullpath'));
    cd(main_path);
    addpath(genpath(fullfile(main_path, 'Algorithms')));
    addpath(genpath(fullfile(main_path, 'Public')));
    addpath(genpath(fullfile(main_path, 'Benchmarks')));
    addpath(genpath(fullfile(main_path, 'Accessories')));
    addpath(genpath(fullfile(main_path, 'GUI')));
end