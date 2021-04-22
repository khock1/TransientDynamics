% Script to run stability analysis/hysteresis diagram for the coral model

% This code supplements the article "Transient dynamics mask the resilience
% of coral reefs" by Hock et al.
% Author & copyright: Karlo Hock, University of Queensland. 2021
% Data should be retrieved from https://zenodo.org/record/4708653

figure;hold on;

a = 0.36;
s = 0.96;
n = 0.51;
g = 0.4;
b = 0.19;
m = 0.04;
h = 0.1;
z = 0.64;
o = 4;
w = 0.4;
f = 1;
p = 1;
k = 0.25;
e = 0.5;
AT= 0:0.05:1;
syms As Ms;
for i=1:length(AT)
   
    e=AT(i);
    
    dA=a*(-(k*p*(e + As*f*w)*(As + Ms - 1))/(a + n + Ms*s + k*p*(e + As*f*w))) + g*As*(1-(-(k*p*(e + As*f*w)*(As + Ms - 1))/(a + n + Ms*s + k*p*(e + As*f*w)))-As-Ms) - b*As*Ms - m*As;
    dM=s*Ms*(1-(-(k*p*(e + As*f*w)*(As + Ms - 1))/(a + n + Ms*s + k*p*(e + As*f*w)))-As-Ms) + s*(-(k*p*(e + As*f*w)*(As + Ms - 1))/(a + n + Ms*s + k*p*(e + As*f*w)))*Ms + b*As*Ms - h*Ms - (z*Ms*o*As)/(1+o*As);
    J=jacobian([dA, dM ],[As Ms]);
    ss1=vpasolve([0 == a*(-(k*p*(e + As*f*w)*(As + Ms - 1))/(a + n + Ms*s + k*p*(e + As*f*w))) + g*As*(1-(-(k*p*(e + As*f*w)*(As + Ms - 1))/(a + n + Ms*s + k*p*(e + As*f*w)))-As-Ms) - b*As*Ms - m*As,0==s*Ms*(1-(-(k*p*(e + As*f*w)*(As + Ms - 1))/(a + n + Ms*s + k*p*(e + As*f*w)))-As-Ms) + s*(-(k*p*(e + As*f*w)*(As + Ms - 1))/(a + n + Ms*s + k*p*(e + As*f*w)))*Ms + b*As*Ms - h*Ms - (z*Ms*o*As)/(1+o*As)],[As,Ms]);
    for jj=1:length(ss1.As)
        if isreal(ss1.As(jj)) && isreal(ss1.Ms(jj))
            if all([ss1.As(jj) ss1.Ms(jj)]>=0 & [ss1.As(jj) ss1.Ms(jj)]<=1)
                jaceigv=eig(double(subs(J,[As Ms],[ss1.As(jj) ss1.Ms(jj)])));
                if any(jaceigv>=0)
                    scatter(AT(i),ss1.As(jj),100,'d','MarkerEdgeColor','k','LineWidth',2.5);%if unstable [0 .7 .7]
                else
                    scatter(AT(i),ss1.As(jj),100,'o','MarkerEdgeColor','k','LineWidth',2.5);%if stable
                end
            end
        end
    end
end
ylabel('Coral Cover');
xlabel('External Supply');
caxis([0 0.05]);
