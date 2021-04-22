% Script to run system trajectories for figure 4

% This code supplements the article "Transient dynamics mask the resilience
% of coral reefs" by Hock et al.
% Author & copyright: Karlo Hock, University of Queensland. 2021

load('results_final_stoch.mat');
scen=1;
relz=50;
ee=[];
thisparams=results_final_stochastic(scen).deltacoral(relz).params;
thisparams(:,14)=results_final_stochastic(scen).deltacoral(relz).params(:,14)/2;
thisparams=vertcat(thisparams,thisparams);
TP=thisparams;


%Fig4A
et=normrnd(0.75,0.01,[100 1]);
et(et>1)=1; 
et(et<0)=0;
thisparams(:,14)=et;
ee(:,1)=et;

Astates = f_run_systemtraj(100, thisparams, 1);

refst=results_final_stochastic(scen).refstates;  
coraltraj=Astates(:,refst);
trj=[6 17];

figure;hold on
for z=1:length(trj)
    i=trj(z);
    plot(1:100, coraltraj(:,i));
end
axis([0 (100+1) 0 1]);
xlabel('Years');
ylabel('Coral Cover');

%Fig4B
et=ee(:,1);
disrupt=25:10:95;
for i=1:length(disrupt)
    thisdl=5;
    et(disrupt(i):(disrupt(i)+thisdl))=normrnd(0.1,0.01,[thisdl+1 1]);
end
et(et>1)=1; 
et(et<0)=0;
ee(:,3)=et;
thisparams(:,14)=et;

Astates = f_run_systemtraj(100, thisparams, 1);

refst=results_final_stochastic(scen).refstates;  
coraltraj=Astates(:,refst);
trj=[11 17];

figure;hold on
for z=1:length(trj)
    i=trj(z);
    plot(1:100, coraltraj(:,i));
end
axis([0 (100+1) 0 1]);
xlabel('Years');
ylabel('Coral Cover');


%Fig4C
et=ee(:,1);
disrupt=5:10:40;
for i=1:length(disrupt)
    thisdl=5;%randi(5)-1;
    et(disrupt(i):(disrupt(i)+thisdl))=normrnd(0.1,0.01,[thisdl+1 1]);
end
et(et>1)=1; 
et(et<0)=0;
ee(:,2)=et;
thisparams(:,14)=et;

Astates = f_run_systemtraj(100, thisparams, 1);

refst=results_final_stochastic(scen).refstates;  
coraltraj=Astates(:,refst);
trj=[6 17];

figure;hold on
for z=1:length(trj)
    i=trj(z);
    plot(1:100, coraltraj(:,i));
end
axis([0 (100+1) 0 1]);
xlabel('Years');
ylabel('Coral Cover');