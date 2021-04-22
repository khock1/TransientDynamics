% Script to reproduce figure 2 trajectories

% This code supplements the article "Transient dynamics mask the resilience
% of coral reefs" by Hock et al.
% Author & copyright: Karlo Hock, University of Queensland. 2021

[Astates, Mstates, col, eql, endst] = f_run_systemtraj(50, thisparams, 1);

refst=resultsscen_final(scen).refstates;  
coraltraj=Astates(:,refst);
trj=sort([6:5:26 17]);

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
ylabel('Adult Coral');