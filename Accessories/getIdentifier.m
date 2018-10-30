function T = getIdentifier(varargin)
if nargin == 0
    time = datetime('now');
elseif nargin == 1
    time = varargin{1};
else
    error('wrong use of timePrefix');
end
T = char(time);
T([strfind(T, '-'), strfind(T, ':')]) = [];
T = strrep(T, ' ', '-');
T(1: 5) = [];
end