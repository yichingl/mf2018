function c = pov( x)
% function c = pov( x)
%  Compute the covariance of x (biased estimate)
%
% Paris Smaragdis, 14-May-03

x = x - repmat( mean( x, 2), 1, cols( x));
c = x*x' / cols( x);
