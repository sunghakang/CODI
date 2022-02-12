function u = normalize (x)
% This is a helper method to normalize a the values in x to [0,255]
% normalize(x)
%     -inputs:
%      x            : Input matrix of size m, n
%     -output :
%      u            : Output matrix of size m, n, values in [0,255]
minx = min(min(x));
x = x - minx;
maxx = max(max(x));
u = (x./maxx)*255;
end