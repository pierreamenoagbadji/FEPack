function xmod = mymod(x, a, b)
% MYMOD(x, a, b) returns xmod such that xmod belongs to (a, b)
% and x - xmod is an integer.
%
% If a = 0, then MYMOD(x, 0, b) simply corresponds to mod(x, b).
if (nargin < 3)
  b = 1;
end

if (nargin < 2)
  a = 0;
end

xmod = a + (b - a) * mod((x - a) / (b - a), 1);
