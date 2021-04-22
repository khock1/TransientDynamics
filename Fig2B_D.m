% Script to reproduce system trajectories for Figure 2

% This code supplements the article "Transient dynamics mask the resilience
% of coral reefs" by Hock et al.
% Author & copyright: Karlo Hock, University of Queensland. 2021

% Data should be retrieved from https://zenodo.org/record/4708653

%Note: Fig2A is a conceptual sketch, not derived from data

relz=50;
scen=1;
load('results_final_determin.mat');%load deterministic results
thisparams=resultsscen_final(scen).deltacoral(relz).params;

%make sure connectivity is zero to get bistability
et=0;
thisparams(:,14)=et;

%Fig2B
reproduce_fig2;

%Fig2C
disrupt=5;
thisparams(disrupt,7)=0.15;
reproduce_fig2;  

%Fig2D
thisparams=resultsscen_final(scen).deltacoral(relz).params;
enhance=5:50;
thisparams(enhance,9)=0.7;
reproduce_fig2;

