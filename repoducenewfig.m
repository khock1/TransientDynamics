% load('results7scen_1k.mat', 'resultsscen')
% load('scenarios.mat')
% resultsscen_final=resultsscen;
scen=1;
relz=50;
modrange=5:10;
trj=sort([6:5:26 17]);
refst=resultsscen_final(scen).refstates;

thisparams=resultsscen_final(scen).deltacoral(relz).params;
thisparams(:,14)=resultsscen_final(scen).deltacoral(relz).params(:,14)/2;

[Astates, Mstates, col, eql, endst] = f_discrER_reps_4figs_02(50, thisparams, 1);


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
% 
% t=8;
% typtraj=resultsscen(scen).typtraj;
% cc=typtraj(:,2);
% mm=typtraj(:,1);
% changeA=zeros(timesteps,length(cc));
% coral=zeros(timesteps,length(cc));
% macroalgae=zeros(timesteps,length(cc));
% coral(1,:)=cc;
% macroalgae(1,:)=mm;
%     for rf=1:length(cc)
%         [changeM(t,rf),changeA(t,rf)]=coral_discER_test([mm(rf) cc(rf)],thisparams(t,:));
%     end
%     coral(t,:)=coral(t-1,:)+changeA(t,:);
%     macroalgae(t,:)=macroalgae(t-1,:)+changeM(t,:);
%     coral(t,coral(t,:)<0)=0;
% 
% figure; hold on
% plot(cc,zeros(length(cc)),'r');
% plot(cc,cumtrapz(-resultsscen(ssc).deltacoral(realiz).changeA(t,:)),'Color','k');
% %plot(cc,resultsscen(ssc).U(realiz,:)/50,'k')
% 
% %ylim([-80 80]);
% xlabel('Coral Cover');
% ylabel('Quasi-potential');
% ylim([-5 5])