function saveGlobal(varargin)
if nargin == 1
    Global = varargin{1};
    directory = 'Data';
elseif nargin == 2
    Global = varargin{1};
    directory = varargin{2};
else
    error('varargin error');
end
if ~exist(directory, 'dir')
    mkdir(directory);
end
T = Global.identifier;
if isempty(obj.algorithm)
    algoname = "unknown";
else
    algoname = obj.algorithm;
end
filepath = sprintf("%s/[%s]%s_%s_F%d.mat", directory, T, algoname, func2str(obj.suite), obj.funcnum);
pointer = 1;
while exist(filepath, 'file')
    filepath = sprintf("%s/[%s]%s_%s_F%d_%d.mat", directory, T, algoname, func2str(obj.suite), obj.funcnum, pointer);
    pointer = pointer + 1;
end
save(filepath, 'Global');
end