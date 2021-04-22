% Script to run the coral  model

% This code supplements the article "Transient dynamics mask the resilience
% of coral reefs" by Hock et al.
% Author & copyright: Karlo Hock, University of Queensland. 2021

function [Mn, An] = coral_discER_backup(y,pr)
p = pr(1);
a = pr(2);
s = pr(3);
n = pr(4);
g = pr(5);
b = pr(6);
m = pr(7);
h = pr(8);
z = pr(9);
o = pr(10);
w = pr(11);
f = pr(12);
k = pr(13);
e = pr(14);
%R = y(1);
A = y(2);
M = y(1);


R=-(k*p*(e + A*f*w)*(A + M - 1))/(a + n + M*s + k*p*(e + A*f*w));
An = a*R + g*A*(1-R-A-M) - b*A*M - m*A;
Mn = s*M*(1-R-A-M) + s*R*M + b*A*M - h*M - (z*M*o*A)/(1+o*A);

end