function T = getIdentifier(varargin)
if nargin == 0
    time = datetime('now');
elseif nargin == 1
    time = varargin{1};
else
    error('wrong use of timePrefix');
end
MO = twoDigitStr(month(time));
DD = twoDigitStr(day(time));
HH = twoDigitStr(hour(time));
MI = twoDigitStr(minute(time));
SS = twoDigitStr(second(time));
T = [MO, DD, '-', HH, MI, SS];
pause(1);
end

function S = twoDigitStr(N)
if N < 10
    S = ['0', num2str(round(N))];
else
    S = num2str(round(N));
end
end