% Script to query sensitivity analysis for finding transients to determine
% proportion fo trajectories with transients

% This code supplements the article "Transient dynamics mask the resilience
% of coral reefs" by Hock et al.
% Author & copyright: Karlo Hock, University of Queensland. 2021

% find if coral attractor exists
% make sure that the coral cover does not go below 5percent in the window
load('sens_traj.mat');
dofigs=0;
windowl=15;
maxdev=0.05;
mincc=0.05;
trj=sort([6:5:26 17]);
transients=zeros(size(store_traj,2),1);
%figure;hold on;
for mcs=1:size(store_traj,2)
    for z=1:length(trj)
        thistraj=store_traj(mcs).coraltraj(:,trj(z));
        eql=store_traj(mcs).eql;
        for stp=1:(50-windowl)
            if (abs(max(thistraj(stp:(stp+(windowl-1)))-thistraj(stp))))<=maxdev%if it does not go above 5%; use moveavg here?
                if (abs(min(thistraj(stp:(stp+(windowl-1)))-thistraj(stp))))<=maxdev%if it does nto go below 5%
                    if min(thistraj(stp:(stp+(windowl-1))))>=mincc%must not drop below 5% total
                        stbeq=vertcat(eql(stp:(stp+(windowl-1))).stabeq);
                        if any(stbeq(:,2)>mincc)%any stab eq greater than 0.05
                            transients(mcs,1)=1;
                        end
                    end
                end
            end
        end
    end
    if dofigs==1
        figure;hold on
        for z=1:length(trj)
            i=trj(z);
            plot(1:50, store_traj(mcs).coraltraj(:,i));
        end
        for tt=1:50
            unseq=store_traj(mcs).eql(tt).unstabeq;
            seq=store_traj(mcs).eql(tt).stabeq;
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
    end
end

figure;
for i=1:length(transients)
    scs(i,1)=sum(transients(1:i),1)/i;
end
plot(1:length(transients),scs);
axis([0 length(transients) 0 1]);
xlabel('Total simulations');
ylabel('Proportion of simulations with transients');