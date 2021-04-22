% Script to run system trajectories for the coral  model

% This code supplements the article "Transient dynamics mask the resilience
% of coral reefs" by Hock et al.
% Author & copyright: Karlo Hock, University of Queensland. 2021

function [Astates, Mstates, col, eql, endstate] = f_run_systemtraj(timesteps, pr, dostab)
col=[];

p = pr(:,1);
a = pr(:,2);
s = pr(:,3);
n = pr(:,4);
g = pr(:,5);
b = pr(:,6);
m = pr(:,7);
h = pr(:,8);
z = pr(:,9);
o = pr(:,10);
w = pr(:,11);
f = pr(:,12);
k = pr(:,13);
e = pr(:,14);

axlim=0.9;

x=0:0.01:axlim;
y=0:0.01:axlim;

spc=[];
for i=1:size(x,2)
    for j=1:size(y,2)
        if x(i)+y(j)<=axlim%<= 1 for sure
            spc=vertcat(spc,[x(i) y(j)]);
        end
    end
end

endstate=zeros(length(unique(spc(:,1))));


xco=0:0.01:1;
yco=0:0.01:1;
htmAf=zeros(length(xco));

Mstates=zeros(timesteps,size(spc,1));
Astates=zeros(timesteps,size(spc,1));

eql=struct('stabeq',[],'unstabeq',[]);
if dostab==1
    [stabeq, unstabeq] = discER_stab([p,a,s,n,g,b,m,h,z,o,w,f,k,e]);
    saddle=[];
    keepsaddle=[];
    for eq=1:size(unstabeq,1)
        if unstabeq(eq,1)>0 && unstabeq(eq,2)>0
            saddle=unstabeq(eq,:);
            keepsaddle=saddle;
        end
    end
    
    for ts=2:timesteps
        [stabeq, unstabeq] = discER_stab([p(ts),a(ts),s(ts),n(ts),g(ts),b(ts),m(ts),h(ts),z(ts),o(ts),w(ts),f(ts),k(ts),e(ts)]);
        eql(ts).stabeq=stabeq;
        eql(ts).unstabeq=unstabeq;
        saddle=[];
        for eq=1:size(unstabeq,1)
            if unstabeq(eq,1)>0 && unstabeq(eq,2)>0
                saddle=unstabeq(eq,:);
                keepsaddle=saddle;
            end
        end
        if isempty(saddle)
            saddle=keepsaddle;
        end
    end
end

for ii = 1:size(spc,1)
    state(1,1:2)=[spc(ii,1),spc(ii,2)];
    for t=2:timesteps
        
        [thisM,thisA]=coral_discER_backup(state(t-1,1:2),[p(t),a(t),s(t),n(t),g(t),b(t),m(t),h(t),z(t),o(t),w(t),f(t),k(t),e(t)]);
        state(t,1)=state(t-1,1)+thisM;
        state(t,2)=state(t-1,2)+thisA;
        if state(t,1)<0
            state(t,1)=0;
        end
        if state(t,2)<0
            state(t,2)=0;
        end
    end
    Mstates(:,ii)=state(:,1);
    Astates(:,ii)=state(:,2);
    
    
end


for ii = 1:size(spc,1)
    
    endstate(x(1,:)==spc(ii,1),y(1,:)==spc(ii,2))=Astates(end,ii);
end

end
