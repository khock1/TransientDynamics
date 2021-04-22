% Script to run system trajectories for figure 3

% This code supplements the article "Transient dynamics mask the resilience
% of coral reefs" by Hock et al.
% Author & copyright: Karlo Hock, University of Queensland. 2021

%Note: Fig3A is a conceptual sketch, not derived from data

load('results_final_stoch.mat');
scen=1;
relz=50;
thisparams=results_final_stochastic(scen).deltacoral(relz).params;

modrange=6:10;%time steps with modified supply
trj=sort([6:5:26 17]);
refst=results_final_stochastic(scen).refstates;  

thisparams=results_final_stochastic(scen).deltacoral(relz).params;
thisparams(:,14)=results_final_stochastic(scen).deltacoral(relz).params(:,14)/2;
et=thisparams(:,14);
et(et>1)=1; 
et(et<0)=0;

[Astates, ~, ~, eql,~] = f_run_systemtraj(50, thisparams, 1);


coraltraj=Astates(:,refst);

%Fig3B
figure;hold on
for z=1:length(trj)
    i=trj(z);
    plot(1:50, coraltraj(:,i));
end

for tt=1:50
    unseq=eql(tt).unstabeq;
    seq=eql(tt).stabeq;
    for jj=1:size(unseq,1)
        if unseq(jj,2)~=0
            scatter(tt,unseq(jj,2),50,'d','k','LineWidth',1.25);%,'filled'
        end
    end
    for jj=1:size(seq,1)
        scatter(tt,seq(jj,2),50,'o','filled','k');%if stable ,'k','filled'
    end
end

axis([0 (50+1) 0 1]);
xlabel('Years');
ylabel('Coral cover');

%Fig3C, supply failure
et3=et;
et3(modrange,1)=0;
thisparams(:,14)=et3;
[Astates, ~, ~, eql, ~] = f_run_systemtraj(50, thisparams, 1);

coraltraj=Astates(:,refst);

figure;hold on
for z=1:length(trj)
    i=trj(z);
    plot(1:50, coraltraj(:,i));
end

for tt=1:50
    unseq=eql(tt).unstabeq;
    seq=eql(tt).stabeq;
    for jj=1:size(unseq,1)
        if unseq(jj,2)~=0
            scatter(tt,unseq(jj,2),50,'d','k','LineWidth',1.25);%,'filled'
        end
    end
    for jj=1:size(seq,1)
        scatter(tt,seq(jj,2),50,'o','filled','k');%if stable ,'k','filled'
    end
end

axis([0 (50+1) 0 1]);
xlabel('Years');
ylabel('Coral cover');

%Fig3D, enhanced supply
et2=et;
et2(modrange,1)=1;
thisparams(:,14)=et2;
[Astates, ~, ~, eql, ~] = f_run_systemtraj(50, thisparams, 1);

coraltraj=Astates(:,refst);

figure;hold on
for z=1:length(trj)
    i=trj(z);
    plot(1:50, coraltraj(:,i));
end

for tt=1:50
    unseq=eql(tt).unstabeq;
    seq=eql(tt).stabeq;
    for jj=1:size(unseq,1)
        if unseq(jj,2)~=0
            scatter(tt,unseq(jj,2),50,'d','k','LineWidth',1.25);%,'filled'
        end
    end
    for jj=1:size(seq,1)
        scatter(tt,seq(jj,2),50,'o','filled','k');%if stable ,'k','filled'
    end
end

axis([0 (50+1) 0 1]);
xlabel('Years');
ylabel('Coral cover');

