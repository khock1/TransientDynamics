% Script to run sensitivity analysis for finding transients

% This code supplements the article "Transient dynamics mask the resilience
% of coral reefs" by Hock et al.
% Author & copyright: Karlo Hock, University of Queensland. 2021

% Data should be retrieved from https://zenodo.org/record/4708653

runs=2400;
store_traj=struct('run',[],'coraltraj',[]);
load('results_final_stoch.mat');
scen=1;
relz=50;
parfor rns=1:runs

    ee=[];
    thisparams=results_final_stochastic(scen).deltacoral(relz).params;
    thisparams(:,14)=results_final_stochastic(scen).deltacoral(relz).params(:,14)/2;
    thisparams=vertcat(thisparams,thisparams);
    thisparams_store=thisparams;
    
    pparams=thisparams;
    pparams(:,2)=(0.55-0.15).*rand(1,1) + 0.15;
    pparams(:,3)=(1-0.75).*rand(1,1) + 0.75;
    pparams(:,4)=(0.7-0.3).*rand(1,1) + 0.3;
    pparams(:,5)=(0.6-0.2).*rand(1,1) + 0.2;
    pparams(:,6)=(0.4-0.01).*rand(1,1) + 0.01;
    nondstc=find(thisparams_store(:,7)==0.04);
    pparams(nondstc,7)=(0.24-0.01).*rand(1,1) + 0.01;
    nondstm=find(thisparams_store(:,8)==0.1);
    pparams(nondstm,8)=(0.3-0.01).*rand(1,1) + 0.01;
    pparams(:,9)=(0.85-0.45).*rand(1,1) + 0.45;
    pparams(:,10)=(6-2).*rand(1,1) + 2;
    pparams(:,11)=(0.6-0.2).*rand(1,1) + 0.2;
    pparams(:,13)=(0.45-0.05).*rand(1,1) + 0.05;

    thisparams=pparams;
    
    [Astates, ~, ~, eql] = f_run_systemtraj(50, thisparams, 1);
    
    refst=results_final_stochastic(scen).refstates;
    coraltraj=Astates(:,refst);
    
    store_traj(rns).run=rns;
    store_traj(rns).coraltraj=coraltraj;
    store_traj(rns).eql=eql;
end
save('sens_traj.mat','store_traj');
