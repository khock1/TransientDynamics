% Script to run stability analysis for the coral  model as a function

% This code supplements the article "Transient dynamics mask the resilience
% of coral reefs" by Hock et al.
% Author & copyright: Karlo Hock, University of Queensland. 2021

function [stabeq, unstabeq] = discER_stab(pr)

syms As Ms;
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

stabeq=[];
unstabeq=[];
dA=a*(-(k*p*(e + As*f*w)*(As + Ms - 1))/(a + n + Ms*s + k*p*(e + As*f*w))) + g*As*(1-(-(k*p*(e + As*f*w)*(As + Ms - 1))/(a + n + Ms*s + k*p*(e + As*f*w)))-As-Ms) - b*As*Ms - m*As;
dM=s*Ms*(1-(-(k*p*(e + As*f*w)*(As + Ms - 1))/(a + n + Ms*s + k*p*(e + As*f*w)))-As-Ms) + s*(-(k*p*(e + As*f*w)*(As + Ms - 1))/(a + n + Ms*s + k*p*(e + As*f*w)))*Ms + b*As*Ms - h*Ms - (z*Ms*o*As)/(1+o*As);
J=jacobian([dA, dM ],[As Ms]);
ss1=vpasolve([0 == a*(-(k*p*(e + As*f*w)*(As + Ms - 1))/(a + n + Ms*s + k*p*(e + As*f*w))) + g*As*(1-(-(k*p*(e + As*f*w)*(As + Ms - 1))/(a + n + Ms*s + k*p*(e + As*f*w)))-As-Ms) - b*As*Ms - m*As,0==s*Ms*(1-(-(k*p*(e + As*f*w)*(As + Ms - 1))/(a + n + Ms*s + k*p*(e + As*f*w)))-As-Ms) + s*(-(k*p*(e + As*f*w)*(As + Ms - 1))/(a + n + Ms*s + k*p*(e + As*f*w)))*Ms + b*As*Ms - h*Ms - (z*Ms*o*As)/(1+o*As)],[As,Ms]);
for jj=1:length(ss1.As)
    if isreal(ss1.As(jj)) && isreal(ss1.Ms(jj))
        if all([ss1.As(jj) ss1.Ms(jj)]>=0 & [ss1.As(jj) ss1.Ms(jj)]<=1)
            jcbeigv=eig(double(subs(J,[As Ms],[ss1.As(jj) ss1.Ms(jj)])));
            if any(jcbeigv>=0)
                unstabeq=double(vertcat(unstabeq,[ss1.Ms(jj) ss1.As(jj)]));
            else
                stabeq=double(vertcat(stabeq,[ss1.Ms(jj) ss1.As(jj)]));
            end
        end
    end
end



end