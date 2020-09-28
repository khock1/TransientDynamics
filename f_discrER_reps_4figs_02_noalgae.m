function [Astates, Mstates, col, eql, endstate] = f_discrER_reps_4figs_02(timesteps, pr, dostab)
%clearvars -except r reps S es ee uui;
%es=[0];
%es=0:0.05:1;
showfigs=[5 20 30];% F1 15 10 20 25 F2 20 F3 5 20 % 1 4 5 8 9 10 11 13 15 16 17 %1 4 5 9 10 15 16
storestates=struct('e',[],'Mstates',[],'Astates',[],'stabeq',[],'unstabeq',[]);
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

% % parameter values ---------
% a = 0.36;%0.35
% s = 0.96;%0.8 0.96
% n = 0.51;%0.45 0.51
% g = 0.4;%0.2 0.4
% b = 0.19;%0.1
% m = 0.04    ;%0.05 078 0.04!
% h = 0.1;%0.2 0.01 
% z = 0.64;%0.64 0.28
% o = 4;%4
% w = 0.4;%0.4
% f = 1;%1
% p = 1;%1
% k = 0.25;%0.1 %0.25!
% e = es(rep);

%stasis modified
% a = 0.36;
% s = 0.6;%8
% n = 0.51;
% g = 0.2;%do not overgrow
% b = 0.19;%1
% m = 0.04;
% h = 0.2;%algae die
% z = 0.05;%not grazed
% o = 4;
% w = 0.4;
% f = 1;
% p = 1;
% k = 0.25;


%stasis
% a = 0.35;
% s = 0.6;%8
% n = 0.45;
% g = 0.2;
% b = 0.1;%1
% m = 0.05;
% h = 0.2;
% z = 0.05;
% o = 4;
% w = 0.7;
% f = 1;
% p = 1;
% k = 0.1;

% a = 0.35;%0.35
% s = 0.8;%0.8
% n = 0.45;%0.45
% g = 0.2;%0.2
% b = 0.1;%0.1
% m = 0.05;%0.05
% h = 0.2;%0.2
% z = 0.4;%0.05
% o = 4;%4
% w = 0.4;%0.4
% f = 1;%1
% p = 1;%1
% k = 0.25;%0.1
% e = es(rep);


% a = 0.35;
% s = 0.8;%8
% n = 0.45;
% g = 0.2;
% b = 0.1;%1
% m = 0.03;
% h = 0.2;
% z = 0.3;
% o = 3;
% w = 0.7;
% f = 1;
% p = 1;
% k = 0.1;
% e = 0.5;


%stasis
% a = 0.35;
% s = 0.6;%8
% n = 0.45;
% g = 0.2;
% b = 0.1;%1
% m = 0.05;
% h = 0.2;
% z = 0.05;
% o = 4;
% w = 0.7;
% f = 1;
% p = 1;
% k = 0.1;
% e = es(rep);

%stasis3 %%% this works really well
% a = 0.35;
% s = 0.96;%8
% n = 0.45;
% g = 0.2;
% b = 0.1;%1
% m = 0.05;
% h = 0.4;
% z = 0.05;
% o = 4;
% w = 0.7;
% f = 1;
% p = 1;
% k = 0.1;
% e = es(rep);

%timesteps=100;
axlim=0.9;

% cyctimes=find(genseqnc('p', timesteps, 7));%find(genseqnc('n', 50, [1 0.1429]));%[find(genseqnc('p', 20, 7)) find(genseqnc('n', 30, [1 0.1429]))];%find(genseqnc('p', timesteps, 7));
% %cyctimes(cyctimes > (timesteps-5))=[];
% %cyctimes=unique([cyctimes 15 17]);
% %cyctimes(cyctimes > 25)=[];
% blchtimes=find(genseqnc('p', timesteps, 7));%find(genseqnc('p', timesteps, 4));
% %blchtimes(blchtimes > (timesteps-5))=[];
% %blchtimes=unique([blchtimes 5 7]);
% graztimes=0;%1:50;
% grazrecovery=210;
% conntimes=1:timesteps;%find(genseqnc('p', timesteps, 2));
% morttimes=0;%1:timesteps;
% hmorttimes=0;%1:timesteps;%this
% cgrowvartimes=0;%1:timesteps;
% spcAgrowtimes=0;%1:timesteps;%this

x=0:0.01:axlim;
y=0:0.01:axlim;
%x=0.25;y=0.25;
spc=[];
for i=1:size(x,2)
    for j=1:size(y,2)
        if x(i)+y(j)<=axlim%<= 1 for sure
            spc=vertcat(spc,[x(i) y(j)]);
        end
    end
end
spc(:,1)=0;
htmA=zeros(length(unique(spc(:,1))));
htmAF=zeros(length(unique(spc(:,1))));
htmM=zeros(length(unique(spc(:,1))));
htmD=zeros(length(unique(spc(:,1))));
htmAP=zeros(length(unique(spc(:,1))));
htmMP=zeros(length(unique(spc(:,1))));
htmBA=zeros(length(unique(spc(:,1))));
htmIGR=zeros(length(unique(spc(:,1))));
htmIGRRT=zeros(length(unique(spc(:,1))));
htmPS=zeros(length(unique(spc(:,1))));
htmTS=zeros(length(unique(spc(:,1))));
htmSLD=zeros(length(unique(spc(:,1))));
htmDFS=zeros(length(unique(spc(:,1))));
endstate=zeros(length(unique(spc(:,1))));
PRM=zeros(length(unique(spc(:,1))),length(e));
CC=zeros(length(unique(spc(:,1))),length(e));
SPD=zeros(length(unique(spc(:,1))),length(e));
%set(gca,'Color','k');
xco=0:0.01:1;
yco=0:0.01:1;
htmAf=zeros(length(xco));
%htmAf=zeros(length(x));
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
    
    % saddle=[0.476,0.169];%saddle if not given
    % if isempty(saddle)
    % 	saddle=[0.424582294814088,0.154044884856278];
    %     %saddle=keepsaddle;
    % end
    
    %[p,a,s,n,g,b,m,h,z,o,w,f,k,e]=dodist(p,a,s,n,g,b,m,h,z,o,w,f,k,e,timesteps,cyctimes,blchtimes,graztimes,conntimes,grazrecovery,morttimes,hmorttimes,cgrowvartimes,spcAgrowtimes);
    
    %if ismember(20,showfigs)
    
    
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
%         if t<11
%             z=0.64;
%         else
%             z=0.05;
%         end
%         if t<11
%             e=1;
%         else
%             e=0.0;
%         end

%         if t==11
%             state((t-1),2)=state((t-1),2)*0.3;
%         end
%         if ismember(t,[8,15])
%             m=0.5;
%             h=0.2;
%         else
%             m = 0.04;
%             h=0.01;
%         end
%         if t<8
%             o=4;
%             g=0.2;
%         else
%             o = 0.1;
%             g=0.1;
%         end
%         %Flood
%         if randi(4)==4
%             s = 0.9;%8
%             n = 0.7;
%             f = 0.5;
%             k = 0.05;
%             g = 0.1;
%             b = 0.2;%1
%         else
%             s = 0.8;%8
%             n = 0.45;
%             f = 1;
%             k = 0.1;
%             g = 0.2;
%             b = 0.1;%1
%         end

          %Recruitment failure
%         if t>10 %mod(t,4)==0 randi(4)==4
%             e=0.05;
%         else
%             e = 0.5;
%         end

        
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
%     v(ii,1)=abs(state(end,2)-state(1,2));
%     v(ii,2)=abs(state(end,1)-state(1,1));
    if ismember(2,showfigs)
        htmA(x(1,:)==spc(ii,1),y(1,:)==spc(ii,2))=abs(state(end,2)-state(1,2));

    end

    if ismember(3,showfigs)
        htmM(x(1,:)==spc(ii,1),y(1,:)==spc(ii,2))=abs(state(end,1)-state(1,1));
    end
    
    if ismember(5,showfigs)
        cg=state(end,2)-state(1,2);
        if cg<0
            cg=0;
        end
        htmAP(x(1,:)==spc(ii,1),y(1,:)==spc(ii,2))=cg;
        htmAF(x(1,:)==spc(ii,1),y(1,:)==spc(ii,2))=state(end,2);
    end
    
    if ismember(6,showfigs)
        ag=state(end,1)-state(1,1);
        if ag<0
            ag=0;
        end
        htmMP(x(1,:)==spc(ii,1),y(1,:)==spc(ii,2))=ag;
    end
    
    if ismember(4,showfigs)
        htmD(x(1,:)==spc(ii,1),y(1,:)==spc(ii,2))=sqrt((state(end,1)-state(1,1))^2+(state(end,2)-state(1,2))^2);
    end
    if ismember(17,showfigs)
        htmDFS(x(1,:)==spc(ii,1),y(1,:)==spc(ii,2))=abs((state(end,1)-state(1,1)))/abs((state(end,2)-state(1,2)));
    end
    %find how many times tthe state was between within certain x,y range

    %htmAf(find(x(1,:)==spc(ii,1)),find(y(1,:)==spc(ii,2)))=nnz(find(state(:,1)>=spc(ii,1) && state(:,1)<spc(ii+1,1)) && find(state(:,2)>=spc(ii,2) && state(:,2)<spc(ii+1,2)));
    %xco=x;yco=y;
    
    if ismember(8,showfigs)
        for xl=2:length(xco)
            lx=xco(xl-1);ux=xco(xl);
            for yl=2:length(yco)
                ly=yco(yl-1);uy=yco(yl);
                htmAf(xco==xco(xl),yco==yco(yl))=htmAf(xco==xco(xl),yco==yco(yl))+sum(state(:,1)>=lx & state(:,1)<ux & state(:,2)>=ly & state(:,2)<uy);
            end
        end
    end
    
    if ismember(7,showfigs)
        %not really basins of attraction, as it does not take into account equilibrium value
        %more like expectation that state will be above value
        if state(end,2)<=0.1
            boa=0;
        else
            boa=1;
        end
        htmBA(x(1,:)==spc(ii,1),y(1,:)==spc(ii,2))=boa;
    end
    
    if ismember(9,showfigs)
        igr=state(2,2)-state(1,2);
        htmIGR(x(1,:)==spc(ii,1),y(1,:)==spc(ii,2))=igr;
    end
    
    if ismember(18,showfigs)
        igrrat=state(2,2)/state(1,2);
        htmIGRRT(x(1,:)==spc(ii,1),y(1,:)==spc(ii,2))=igrrat;
    end
    
    if ismember(12,showfigs)
        %how close it passes to saddle point
        thisdist=0;
        for dte=1:size(Astates,1)
            if sqrt((state(dte,1)-saddle(1))^2+(state(dte,2)-saddle(2))^2)<=0.2
                thisdist=1;
                break;
            end
        end
        if thisdist==1
            td=1;
        else
            td=0;
        end
        htmPS(x(1,:)==spc(ii,1),y(1,:)==spc(ii,2))=td;
    end
    
    if ismember(13,showfigs)
        %how long does it spend around saddle point?
        thistime=0;
        for dte=1:size(Astates,1)
            if sqrt((state(dte,1)-saddle(1))^2+(state(dte,2)-saddle(2))^2)<=0.1
                thistime=thistime+1;
            end
        end
        htmTS(x(1,:)==spc(ii,1),y(1,:)==spc(ii,2))=thistime;
    end
    
    if ismember(14,showfigs)
        %does it ever expereicne slowdown beyond some threshold
        sldany=0;
        for st=1:(size(Astates,1)-1)
            if sqrt(state(st+1,1)-(state(st,1))^2+(state(st+1,2)-state(st,2))^2)<=0.05
                sldany=1;
                break;
            end
        end
        if sldany==1
            sld=1;
        else
            sld=0;
        end
        htmSLD(x(1,:)==spc(ii,1),y(1,:)==spc(ii,2))=sld;
    end
    
end

% storestates(rep).e=es(rep);
% storestates(rep).Mstates=Mstates;
% storestates(rep).Astates=Astates;
% storestates(rep).stabeq=stabeq;
% storestates(rep).unstabeq=unstabeq;

%scatter(0.476,0.169);%saddle

%phase plane trajectories
if ismember(1,showfigs)
    figure;hold;
    for iii = 1:size(spc,1)
        if mod(iii,5)==0
            plot(Mstates(:,iii),Astates(:,iii),'Color','k');
%             plot(Mstates(end,iii),Astates(end,iii),'ro')
        end
    end
%     for jj=1:size(unstabeq,1)
%         scatter(unstabeq(jj,1),unstabeq(jj,2),140,'d','m','filled');
%     end
%     for jj=1:size(stabeq,1)
%         scatter(stabeq(jj,1),stabeq(jj,2),140,'d','g','filled');%if stable
%     end
    ylabel('Adult Coral');
    xlabel('Macroalgae');
    axis([0 axlim 0 axlim]);
end

% syms As Ms
% R=-(k*p*(e + As*f*w)*(As + Ms - 1))/(a + n + Ms*s + k*p*(e + As*f*w));
% dy1 = a*R + g*As*(1-R-As-Ms) - b*As*Ms - m*As;
% dy2 = s*Ms*(1-R-As-Ms) + s*R*Ms + b*As*Ms - h*Ms - (z*Ms*o*As)/(1+o*As);
% S1=solve(dy1==0, As);
% S2=solve(dy2==0, As);
% % plh=fplot(S1(1),[0 1]);set(plh,'color','m','linewidth',2);
% % plh=fplot(S2(1),[0 1]);set(plh,'color','green','linewidth',2);
% plh=fplot(S1(2),[0 1]);set(plh,'color','b','linewidth',2);
% plh=fplot(S2(2),[0 1]);set(plh,'color','green','linewidth',2);
% axis([0 axlim 0 axlim]);

%absolute change in macroalgae start to end
if ismember(2,showfigs)
    figure;
    imagesc(x,y,htmM);
    view(-90, 90);
    colormap jet;
    %colorbar;
    xlabel('Adult Coral');
    ylabel('Macroalgae');
end

%absolute change in coral start to end
if ismember(3,showfigs)
    figure;
    imagesc(x,y,htmA)
    view(-90, 90);
    colormap jet;
    %colorbar
    xlabel('Adult Coral');
    ylabel('Macroalgae');
end

%absolute change in system state start to end
if ismember(4,showfigs)
    figure;
    imagesc(x,y,htmD,[0 0.5])
    axis([0 axlim 0 axlim]);
    view(-90, 90);
    %colormap jet;
    colorbar
    xlabel('Adult Coral');
    ylabel('Macroalgae');
end

%positive change in coral state start to end
if ismember(5,showfigs)
%     figure;
%     imagesc(x,y,htmAP,[0 0.5])
%     axis([0 axlim 0 axlim]);
%     view(-90, 90);
%     colormap jet;
%     %colorbar
%     xlabel('Adult Coral');
%     ylabel('Macroalgae');
end

%positive change in macroalgal state start to end
if ismember(6,showfigs)
    figure;
    imagesc(x,y,htmMP)
    view(-90, 90);
    colormap jet;
    %colorbar
    xlabel('Adult Coral');
    ylabel('Macroalgae');
end

%whether coral finishes above some exptected value
if ismember(7,showfigs)
    figure;
    imagesc(x,y,htmBA)
    axis([0 axlim 0 axlim]);
    view(-90, 90);
    colormap jet;
    xlabel('Adult Coral');
    ylabel('Macroalgae');
end

%how often a system is found in a given state
if ismember(8,showfigs)
    figure;
    htmAf=htmAf./max(max(htmAf));
    imagesc(x,y,htmAf,[0 0.005])
    axis([0 axlim 0 axlim]);
    view(-90, 90);
    colormap jet;
    %colorbar
    xlabel('Adult Coral');
    ylabel('Macroalgae');
end

%system change after one timestep
if ismember(9,showfigs)
    figure;
    imagesc(x,y,htmIGR,[0 0.05])
    axis([0 axlim 0 axlim]);
    view(-90, 90);
    colormap jet;
    colorbar
    xlabel('Adult Coral');
    ylabel('Macroalgae');
end

stt1=find(Astates(1,:)==0.24 & Mstates(1,:)==0.24);
stt2=find(Astates(1,:)==0.1 & Mstates(1,:)==0.3);
stt3=find(Astates(1,:)==0.04 & Mstates(1,:)==0.3);
stt4=find(Astates(1,:)==0.1 & Mstates(1,:)==0.1);
stt5=find(Astates(1,:)==0.02 & Mstates(1,:)==0.02);
stt6=find(Astates(1,:)==0.04 & Mstates(1,:)==0.02);
stt7=find(Astates(1,:)==0.16 & Mstates(1,:)==0.4);
stt8=find(Astates(1,:)==0.18 & Mstates(1,:)==0.27);% 14 32
stt9=find(Astates(1,:)==0.05 & Mstates(1,:)==0.5);
stt10=find(Astates(1,:)==0.1 & Mstates(1,:)==0.4);
stt11=find(Astates(1,:)==0.05 & Mstates(1,:)==0.5);
stt12=602;%436
stts=[stt1 stt2 stt3 stt4 stt5 stt6 stt7 stt8 stt9 stt10 stt11 stt12];

% dd1=find(Astates(1,:)==0.1 & Mstates(1,:)==0.1);
% dd2=find(Astates(1,:)==0.1 & Mstates(1,:)==0.4);
% dd3=find(Astates(1,:)==0.1 & Mstates(1,:)==0.7);
% dd4=find(Astates(1,:)==0.2 & Mstates(1,:)==0.2);
% dd5=find(Astates(1,:)==0.2 & Mstates(1,:)==0.5);
% dd6=find(Astates(1,:)==0.3 & Mstates(1,:)==0.3);
% dd1=find(Astates(1,:)==0.1 & Mstates(1,:)==0.2);
% dd2=find(Astates(1,:)==0.12 & Mstates(1,:)==0.4);%0.12 0.4
% dd3=find(Astates(1,:)==0.1 & Mstates(1,:)==0.5);%0.13 0.5
% dd4=find(Astates(1,:)==0.25 & Mstates(1,:)==0.25);
% dd5=find(Astates(1,:)==0.2 & Mstates(1,:)==0.4);
% dd6=find(Astates(1,:)==0.05 & Mstates(1,:)==0.65);
% dd7=find(Astates(1,:)==0.05 & Mstates(1,:)==0.5);
% dds=[dd1 dd2 dd3 dd4 dd5 dd6];
dds=[4005;3717;3413;3117;2825;2557;2267];

%plot system change trajectories
if ismember(10,showfigs)
    figure;hold
%     plot(Mstates(:,stt1),Astates(:,stt1))
%     plot(Mstates(:,stt2),Astates(:,stt2))
%     plot(Mstates(:,stt3),Astates(:,stt3))
%     plot(Mstates(:,stt4),Astates(:,stt4))
%     plot(Mstates(:,stt5),Astates(:,stt5))
%     plot(Mstates(:,stt6),Astates(:,stt6))
%     plot(Mstates(:,stt7),Astates(:,stt7))
%     plot(Mstates(:,stt9),Astates(:,stt9))
%     plot(Mstates(:,stt10),Astates(:,stt10))
%     plot(Mstates(:,stt11),Astates(:,stt11))
%     plot(Mstates(:,stt12),Astates(:,stt12))
    %plot(Mstates(:,stt8),Astates(:,stt8))
    %scatter(saddle(1),saddle(2),140,'d','m');
    
    %plot(Mstates(:,dd1),Astates(:,dd1))
    plot(Mstates(:,dd2),Astates(:,dd2))
    plot(Mstates(:,dd3),Astates(:,dd3))
    plot(Mstates(:,dd4),Astates(:,dd4))
    plot(Mstates(:,dd5),Astates(:,dd5))
    plot(Mstates(:,dd6),Astates(:,dd6))
    %plot(Mstates(:,dd7),Astates(:,dd7))
    
    for jj=1:size(unstabeq,1)
        scatter(unstabeq(jj,1),unstabeq(jj,2),150,'+','k','LineWidth',2.5);
    end
    for jj=1:size(stabeq,1)
        scatter(stabeq(jj,1),stabeq(jj,2),150,'x','k','LineWidth',2.5);%if stable
    end
    axis([0 axlim 0 axlim]);
    ylabel('Adult Coral');
    xlabel('Macroalgae');
end

%plot coral recovery trajectories
if ismember(11,showfigs)
    figure;hold on;
%     plot(1:length(Astates(:,stt1)),Astates(:,stt1))
%      plot(1:length(Astates(:,stt2)),Astates(:,stt2))
%     plot(1:length(Astates(:,stt3)),Astates(:,stt3))
%     plot(1:length(Astates(:,stt4)),Astates(:,stt4))
%     plot(1:length(Astates(:,stt5)),Astates(:,stt5))
%     plot(1:length(Astates(:,stt6)),Astates(:,stt6))
%     plot(1:length(Astates(:,stt7)),Astates(:,stt7))
%     plot(1:length(Astates(:,stt9)),Astates(:,stt9))
%     plot(1:length(Astates(:,stt10)),Astates(:,stt10))
%     plot(1:length(Astates(:,stt11)),Astates(:,stt11))
%     plot(1:length(Astates(:,stt12)),Astates(:,stt12))
    %plot(1:length(Astates(:,stt8)),Astates(:,stt8))
    for i=1:length(dds)
        plot(1:length(Astates(:,1)),Astates(:,dds(i)));
    end
    axis([0 (length(Astates(:,stt1))+1) 0 1]);
    xlabel('Years');
    ylabel('Adult Coral');
end

%plot macroalgal recovery trajectories
if ismember(19,showfigs)
    figure;hold;
    plot(1:length(Mstates(:,stt1)),Mstates(:,stt1))
    plot(1:length(Mstates(:,stt2)),Mstates(:,stt2))
    plot(1:length(Mstates(:,stt3)),Mstates(:,stt3))
    plot(1:length(Mstates(:,stt4)),Mstates(:,stt4))
    plot(1:length(Mstates(:,stt5)),Mstates(:,stt5))
    plot(1:length(Mstates(:,stt6)),Mstates(:,stt6))
    plot(1:length(Mstates(:,stt7)),Mstates(:,stt7))
    plot(1:length(Mstates(:,stt9)),Mstates(:,stt9))
    plot(1:length(Mstates(:,stt10)),Mstates(:,stt10))
    plot(1:length(Mstates(:,stt11)),Mstates(:,stt11))
    plot(1:length(Mstates(:,stt12)),Mstates(:,stt12))
    %plot(1:length(Mstates(:,stt8)),Mstates(:,stt8))
    axis([0 (length(Mstates(:,stt1))+1) 0 1]);
    xlabel('Years');
    ylabel('Macroalgae');
end



%whether system trajectory passes wihtin some distance of saddle - bugged?
if ismember(12,showfigs)
    figure;
    imagesc(x,y,htmPS)
    axis([0 axlim 0 axlim]);
    view(-90, 90);
    colormap jet;
    xlabel('Adult Coral');
    ylabel('Macroalgae');
end

%how long the system trajectory spends within some distance of the saddle point
if ismember(13,showfigs)
    figure;
    imagesc(x,y,htmTS,[0 20])
    axis([0 axlim 0 axlim]);
    view(-90, 90);
    colormap jet;
    %colorbar
    xlabel('Adult Coral');
    ylabel('Macroalgae');
end

%does trajectory ever expereicnes slowdown beyond some threshold - bugged?
if ismember(14,showfigs)
    figure;
    imagesc(x,y,htmSLD)
    axis([0 axlim 0 axlim]);
    view(-90, 90);
    colormap jet;
    xlabel('Adult Coral');
    ylabel('Macroalgae');
end

%nullclines and equilibria
if ismember(15,showfigs)
%     syms As Ms
%     R=-(k*p*(e + As*f*w)*(As + Ms - 1))/(a + n + Ms*s + k*p*(e + As*f*w));
%     dy1 = a*R + g*As*(1-R-As-Ms) - b*As*Ms - m*As;
%     dy2 = s*Ms*(1-R-As-Ms) + s*R*Ms + b*As*Ms - h*Ms - (z*Ms*o*As)/(1+o*As);
%     S1=solve(dy1==0, As);
%     S2=solve(dy2==0, As);
    syms As Ms
    R=-(k(1)*p(1)*(e(1) + As*f(1)*w(1))*(As + Ms - 1))/(a(1) + n(1) + Ms*s(1) + k(1)*p(1)*(e(1) + As*f(1)*w(1)));
    dy1 = a(1)*R + g(1)*As*(1-R-As-Ms) - b(1)*As*Ms - m(1)*As;
    dy2 = s(1)*Ms*(1-R-As-Ms) + s(1)*R*Ms + b(1)*As*Ms - h(1)*Ms - (z(1)*Ms*o(1)*As)/(1+o(1)*As);
    S1=solve(dy1==0, As);
    S2=solve(dy2==0, As);
    figure; hold;
    plh=fplot(S1(1),[0 1]);set(plh,'color','b','linewidth',2);%plh=plot(0:0.1:1,zeros(1,11));set(plh,'color','b','linewidth',2);plh=fplot(S1(1),[0 1]);set(plh,'color','b','linewidth',2);
    %plh=fplot(S2(1),[0 1]);set(plh,'color','green','linewidth',2);
    plh=plot(zeros(1,11),0:0.1:1);set(plh,'color',[0.9290, 0.6940, 0.1250],'linewidth',2);%plh=plot(zeros(1,11),0:0.1:1);set(plh,'color','green','linewidth',2);plh=fplot(S2(1),[0 1]);set(plh,'color','green','linewidth',2);
    plh=fplot(S1(2),[0 1]);set(plh,'color','b','linewidth',2);
    plh=fplot(S2(2),[0 1]);set(plh,'color',[0.9290, 0.6940, 0.1250],'linewidth',2);
    for jj=1:size(unstabeq,1)
        scatter(unstabeq(jj,1),unstabeq(jj,2),150,'+','k','LineWidth',2.5);
    end
    for jj=1:size(stabeq,1)
        scatter(stabeq(jj,1),stabeq(jj,2),150,'x','k','LineWidth',2.5);%if stable
    end
    %scatter(Mstates(1,dd1),Astates(1,dd1),50,'k','filled');
%     scatter(Mstates(1,dd2),Astates(1,dd2),50,'k','filled');
%     scatter(Mstates(1,dd3),Astates(1,dd3),50,'k','filled');
%     scatter(Mstates(1,dd4),Astates(1,dd4),50,'k','filled');
%     scatter(Mstates(1,dd5),Astates(1,dd5),50,'k','filled');
%     scatter(Mstates(1,dd6),Astates(1,dd6),50,'k','filled');
    %scatter(Mstates(1,dd7),Astates(1,dd7),50,'k','filled');
    ylabel('Adult Coral');
    xlabel('Macroalgae');
    axis([0 axlim 0 axlim]);
end

%just equilibria for overlaying on figures
if ismember(16,showfigs)
    figure;hold;
    for jj=1:size(unstabeq,1)
        scatter(unstabeq(jj,1),unstabeq(jj,2),140,'d','m','filled');
    end
    for jj=1:size(stabeq,1)
        scatter(stabeq(jj,1),stabeq(jj,2),140,'d','g','filled');%if stable
    end
%     if isempty(unstabeq)
%         scatter(saddle(1),saddle(2),140,'d','m');
%     end
    ylabel('Adult Coral');
    xlabel('Macroalgae');
    axis([0 axlim 0 axlim]);
end

%ratio of total algal change over total coral change - positive means algae changed nmuch more
if ismember(17,showfigs)
    figure;
    imagesc(x,y,htmDFS,[0 10])
    axis([0 axlim 0 axlim]);
    view(-90, 90);
    colormap jet;
    %colorbar
    xlabel('Adult Coral');
    ylabel('Macroalgae');
end

%ratio of 1-yr algal change over 1-yr coral change - positive means algae changed much more
if ismember(18,showfigs)
    figure;
    imagesc(x,y,htmIGRRT)
    axis([0 axlim 0 axlim]);
    view(-90, 90);
    colormap jet;
    %colorbar
    xlabel('Adult Coral');
    ylabel('Macroalgae');
end

%recovery trajectories and equilibria over time
if ismember(20,showfigs)
%     figure;hold on;
% %     plot(1:length(Astates(:,stt1)),Astates(:,stt1))
% %     plot(1:length(Astates(:,stt2)),Astates(:,stt2))
% %     plot(1:length(Astates(:,stt3)),Astates(:,stt3))
% %     plot(1:length(Astates(:,stt4)),Astates(:,stt4))
% %     plot(1:length(Astates(:,stt5)),Astates(:,stt5))
% %     plot(1:length(Astates(:,stt6)),Astates(:,stt6))
% %     plot(1:length(Astates(:,stt7)),Astates(:,stt7))
% %     plot(1:length(Astates(:,stt9)),Astates(:,stt9))
% %     plot(1:length(Astates(:,stt10)),Astates(:,stt10))
% %     plot(1:length(Astates(:,stt11)),Astates(:,stt11))
% %     plot(1:length(Astates(:,stt12)),Astates(:,stt12))
%     
%     %plot(1:length(Astates(:,dd1)),Astates(:,dd1))
% %     plot(1:length(Astates(:,dd2)),Astates(:,dd2))
% %     plot(1:length(Astates(:,dd3)),Astates(:,dd3))
% %     plot(1:length(Astates(:,dd4)),Astates(:,dd4))
% %     plot(1:length(Astates(:,dd5)),Astates(:,dd5))
% %     plot(1:length(Astates(:,dd6)),Astates(:,dd6))
%     %plot(1:length(Astates(:,dd7)),Astates(:,dd7))
%     for i=1:length(dds)
%         plot(1:length(Astates(:,1)),Astates(:,dds(i)));
%     end
%     
%     
%     %plot(1:length(Astates(:,stt8)),Astates(:,stt8))
%     ue=[];
%     for tt=1:timesteps
%         unseq=eql(tt).unstabeq;
%         seq=eql(tt).stabeq;
%         for jj=1:size(unseq,1)
%             if unseq(jj,2)~=0
%                 scatter(tt,unseq(jj,2),50,'+','k','LineWidth',1.25);%,'filled'
%                 ue=[ue unseq(jj,2)];
%             end
%         end
%         for jj=1:size(seq,1)
%             scatter(tt,seq(jj,2),50,'x','k','LineWidth',1.25);%if stable ,'k','filled'
%         end
% %         if isempty(unseq)
% %             scatter(saddle(1),saddle(2),140,'d','m');
% %         end
%     end
%     axis([0 (length(Astates(:,stt1))+1) 0 1]);
%     xlabel('Years');
%     ylabel('Adult Coral');
%     legend;
end

%IGR
thisst=stt11;
if ismember(21,showfigs)
    figure;hold;
    for tt=1:(timesteps-1)
        igr(tt)=log((Astates(tt+1,thisst)+0.05)/(Astates(tt,thisst)+0.05))/((tt+1)-tt);
    end
    plot(2:length(Astates(:,thisst)),igr);
    %axis([0 (length(Astates(:,stt1))+1) 0 1]);
    xlabel('Years');
    ylabel('IGR');
end

%lambda
if ismember(22,showfigs)
    figure;hold;
    lambda=zeros(timesteps,1);
    for lm=2:timesteps
        lambda(lm,1)=(sum(Astates(:,thisst)))^(1/(lm));
    end
    plot(1:length(Astates(:,thisst)),lambda);
    xlabel('Years');
    ylabel('Lambda');
end

%macroalgal recovery trajectories
if ismember(23,showfigs)
    figure;hold;
    plot(1:length(Mstates(:,stt1)),Mstates(:,stt1))
    plot(1:length(Mstates(:,stt2)),Mstates(:,stt2))
    plot(1:length(Mstates(:,stt3)),Mstates(:,stt3))
    plot(1:length(Mstates(:,stt4)),Mstates(:,stt4))
    plot(1:length(Mstates(:,stt5)),Mstates(:,stt5))
    plot(1:length(Mstates(:,stt6)),Mstates(:,stt6))
    plot(1:length(Mstates(:,stt7)),Mstates(:,stt7))
    plot(1:length(Mstates(:,stt9)),Mstates(:,stt9))
    plot(1:length(Mstates(:,stt10)),Mstates(:,stt10))
    plot(1:length(Mstates(:,stt11)),Mstates(:,stt11))
    plot(1:length(Mstates(:,stt12)),Mstates(:,stt12))



    %plot(1:length(Mstates(:,stt8)),Mstates(:,stt8))
%     for tt=1:timesteps
%         unseq=eql(tt).unstabeq;
%         seq=eql(tt).stabeq;
%         for jj=1:size(unseq,1)
%             scatter(tt,unseq(jj,2),140,'d','m','filled');
%         end
%         for jj=1:size(seq,1)
%             scatter(tt,seq(jj,2),140,'d','g','filled');%if stable
%         end
% %         if isempty(unseq)
% %             scatter(saddle(1),saddle(2),140,'d','m');
% %         end
%     end
    axis([0 (length(Mstates(:,stt1))+1) 0 1]);    
    ylabel('Adult Coral');
    xlabel('Macroalgae');
    legend;
end

%amoutn fo external suppyl per timestep
if ismember(24,showfigs)
    figure;hold;
    scatter(1:timesteps,e);
end

%multicolour phase plot
if ismember(25,showfigs)
    figure;hold;
    for iii = 1:size(spc,1)
        if mod(iii,10)==0
            thissst=iii;
            xx=transpose(Mstates(1:(end-1),thissst));
            yy=transpose(Astates(1:(end-1),thissst));
            ch=zeros((length(xx)),1);
            for c=1:(length(xx))
                ch(c)=pdist2([Mstates(c,thissst) Astates(c,thissst)],[Mstates(c+1,thissst) Astates(c+1,thissst)],'euclidean');
            end
            zz = zeros(size(xx));
            col=transpose(ch);
            surface([xx;xx],[yy;yy],[zz;zz],[col;col],...
                'facecol','no',...
                'edgecol','interp',...
                'linew',1);
            
        end

        
    end
    colormap jet;
%     scatter(Mstates(1,dd1),Astates(1,dd1),50,'k','filled');
%     scatter(Mstates(1,dd2),Astates(1,dd2),50,'k','filled');
%     scatter(Mstates(1,dd3),Astates(1,dd3),50,'k','filled');
%     scatter(Mstates(1,dd4),Astates(1,dd4),50,'k','filled');
%     scatter(Mstates(1,dd5),Astates(1,dd5),50,'k','filled');
%     scatter(Mstates(1,dd6),Astates(1,dd6),50,'k','filled');
%     scatter(Mstates(1,dd7),Astates(1,dd7),50,'k','filled');
    %%draw specific trajectories
%     for iii = 1:size(spc,1)
%         if ismember(iii,stts)
%             plot(Mstates(:,iii),Astates(:,iii),'Color','k','linew',1);
%         end
%     end

    %%draw nullclines
%     syms As Ms
%     R=-(k(1)*p(1)*(e(1) + As*f(1)*w(1))*(As + Ms - 1))/(a(1) + n(1) + Ms*s(1) + k(1)*p(1)*(e(1) + As*f(1)*w(1)));
%     dy1 = a(1)*R + g(1)*As*(1-R-As-Ms) - b(1)*As*Ms - m(1)*As;
%     dy2 = s(1)*Ms*(1-R-As-Ms) + s(1)*R*Ms + b(1)*As*Ms - h(1)*Ms - (z(1)*Ms*o(1)*As)/(1+o(1)*As);
%     S1=solve(dy1==0, As);
%     S2=solve(dy2==0, As);
%     plh=fplot(S1(1),[0 1]);set(plh,'color','black','linewidth',2);%plh=plot(0:0.1:1,zeros(1,11));set(plh,'color','b','linewidth',2);plh=fplot(S1(1),[0 1]);set(plh,'color','b','linewidth',2);
%     plh=plot(zeros(1,11),0:0.1:1);set(plh,'color','black','linewidth',2);%plh=plot(zeros(1,11),0:0.1:1);set(plh,'color','green','linewidth',2);plh=fplot(S2(1),[0 1]);set(plh,'color','green','linewidth',2);
%     plh=fplot(S1(2),[0 1]);set(plh,'color','black','linewidth',2);
%     plh=fplot(S2(2),[0 1]);set(plh,'color','black','linewidth',2);
%     
%     for jj=1:size(unstabeq,1)
%         scatter(unstabeq(jj,1),unstabeq(jj,2),150,'+','k','LineWidth',2.5);
%     end
%     for jj=1:size(stabeq,1)
%         scatter(stabeq(jj,1),stabeq(jj,2),150,'x','k','LineWidth',2.5);%if stable
%     end
    ylabel('Adult Coral');
    xlabel('Macroalgae');
    axis([0 axlim 0 axlim]);

end

%how often coral above threhsold, i.e., how fast coral reaches threshold
if ismember(26,showfigs)
    thresh=0.8;
    for ii = 1:size(spc,1)
        state(:,1)=Mstates(:,ii);
        state(:,2)=Astates(:,ii);
        aat=nnz(find(state(:,2)>=thresh));
        endstate(x(1,:)==spc(ii,1),y(1,:)==spc(ii,2))=aat;
    end
    
    figure;
    imagesc(x,y,endstate)

    axis([0 axlim 0 axlim]);
    view(-90, 90);
    hold;
    %scatter(Astates(1,dd1),Mstates(1,dd1),50,'w','filled');
    scatter(Astates(1,dd2),Mstates(1,dd2),50,'w','filled');
    scatter(Astates(1,dd3),Mstates(1,dd3),50,'w','filled');
    scatter(Astates(1,dd4),Mstates(1,dd4),50,'w','filled');
    scatter(Astates(1,dd5),Mstates(1,dd5),50,'w','filled');
    scatter(Astates(1,dd6),Mstates(1,dd6),50,'w','filled');
    %scatter(Astates(1,dd7),Mstates(1,dd7),50,'w','filled');
    colormap jet;
    %colorbar
    xlabel('Adult Coral');
    ylabel('Macroalgae');

end


%hysteresis speed diagram; only do this is you have chain repetitions
if ismember(27,showfigs)
    %I need this for hysteresis speed diagram
    %I also need to turn on fig 8 and 9 to get htmAf and htmIGR
    trjXM=zeros(length(htmAf),1);
    trjYC=zeros(length(htmAf),1);
    for i=2:length(htmAf)
        trjYC(i,1)=yco(i);
        thismax=find(htmAf(:,i)==max(htmAf(:,i)));
        trjXM(i,1)=xco(thismax(1));
    end
    fz=find(trjYC==0);
    trjYC(fz)=NaN;
    trjXM(fz)=NaN;
    %figure;
    %plot(trjXM(2:82),trjYC(2:82))
    
    trjZS=zeros(length(htmIGR),1);
    for i=1:length(htmIGR)
        if ~isnan(trjXM(i,1))
            if ~isnan(trjYC(i,1))
                trjZS(i,1)=htmIGR(find(xco==trjXM(i,1)),find(yco==trjYC(i,1)));
            end
        end
    end
    PRM(:,rep)=es(rep);
    CC(:,rep)=trjYC(2:length(unique(spc(:,1))));
    SPD(:,rep)=trjZS;
end

% if rep~=length(e)
%     clearvars -except e showfigs storestates keepsaddle;
% end


%BW phase plot
if ismember(28,showfigs)
    figure;hold;
    for iii = 1:size(spc,1)
        if mod(iii,10)==0
            if abs(0.9-Astates(end,iii))<0.1
                plot(Mstates(:,iii),Astates(:,iii),'Color','k');
            elseif abs(0-Astates(end,iii))<0.1
                plot(Mstates(:,iii),Astates(:,iii),'Color',[0.7 0.97 0.7]);
            else
                plot(Mstates(:,iii),Astates(:,iii),'Color',[0.9 0.9 0.9]);
            end
        end

    end
    for jj=1:size(unstabeq,1)
        scatter(unstabeq(jj,1),unstabeq(jj,2),150,'+','r','LineWidth',2.5);
    end
    for jj=1:size(stabeq,1)
        scatter(stabeq(jj,1),stabeq(jj,2),150,'x','r','LineWidth',2.5);%if stable
    end
    ylabel('Adult Coral');
    xlabel('Macroalgae');
    axis([0 axlim 0 axlim]);

end
%how often coral above threhsold, i.e., how fast coral reaches threshold;
%for batch wedge plot
if ismember(30,showfigs)
    %thresh=0.5;
    for ii = 1:size(spc,1)
%         state(:,1)=Mstates(:,ii);
%         state(:,2)=Astates(:,ii);
        %aat=nnz(find(state(:,2)>=thresh));
        endstate(x(1,:)==spc(ii,1),y(1,:)==spc(ii,2))=Astates(end,ii);
    end
%     
%     figure;
%     imagesc(x,y,endstate)
% 
%     axis([0 axlim 0 axlim]);
%     view(-90, 90);
%     hold;
%     %scatter(Astates(1,dd1),Mstates(1,dd1),50,'w','filled');
% %     scatter(Astates(1,dd2),Mstates(1,dd2),50,'w','filled');
% %     scatter(Astates(1,dd3),Mstates(1,dd3),50,'w','filled');
% %     scatter(Astates(1,dd4),Mstates(1,dd4),50,'w','filled');
% %     scatter(Astates(1,dd5),Mstates(1,dd5),50,'w','filled');
% %     scatter(Astates(1,dd6),Mstates(1,dd6),50,'w','filled');
%     %scatter(Astates(1,dd7),Mstates(1,dd7),50,'w','filled');
%     colormap jet;
%     %colorbar
%     xlabel('Adult Coral');
%     ylabel('Macroalgae');

end


if ismember(27,showfigs)
    figure;
    sf = pcolor(PRM,CC,abs(SPD));
    sf.EdgeColor = 'none';
    sf.FaceColor = 'interp';
end
end
% htmAPdif=htmAPold-htmAP;
% figure;
% imagesc(x,y,htmAPdif,[0 0.16])
% view(-90, 90);