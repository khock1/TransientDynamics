% Script to create phase plots (binary and heatmap) for the coral  model

% This code supplements the article "Transient dynamics mask the resilience
% of coral reefs" by Hock et al.
% Author & copyright: Karlo Hock, University of Queensland. 2021


es=[0 0.8];
storestates=struct('e',[],'Mstates',[],'Astates',[],'stabeq',[],'unstabeq',[]);
for rep=1:length(es)
    % parameter values ---------
    a = 0.35;
    s = 0.95;
    n = 0.5;
    g = 0.4;
    b = 0.2;
    m = 0.05;
    h = 0.1;
    z = 0.65;
    o = 4;
    w = 0.4;
    f = 1;
    p = 1;
    k = 0.25;
    e = es(rep);
    
    %disturbance parameters; not used right now
    timesteps=100;
    axlim=0.9;
    cyctimes=0;
    blchtimes=0;
    graztimes=0;
    grazrecovery=210;
    conntimes=0;
    morttimes=0;
    
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


    Mstates=zeros(timesteps,size(spc,1));
    Astates=zeros(timesteps,size(spc,1));
    
    %calculate equilibria
    [stabeq, unstabeq] = discER_stab([p,a,s,n,g,b,m,h,z,o,w,f,k,e]);
    saddle=[];
    keepsaddle=[];
    for eq=1:size(unstabeq,1)
        if unstabeq(eq,1)>0 && unstabeq(eq,2)>0
            saddle=unstabeq(eq,:);
            keepsaddle=saddle;
        end
    end
    
    [p,a,s,n,g,b,m,h,z,o,w,f,k,e]=dodist(p,a,s,n,g,b,m,h,z,o,w,f,k,e,timesteps,cyctimes,blchtimes,graztimes,conntimes,grazrecovery,morttimes,0, 0, 0);
    
    %obtain state trajectories
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
    
    %two-colour phase plot
    figure;hold;
    for iii = 1:size(spc,1)
        if mod(iii,10)==0
            if abs(0.9-Astates(end,iii))<0.1
                plot(Mstates(:,iii),Astates(:,iii),'Color','r');
            elseif abs(0-Astates(end,iii))<0.1
                plot(Mstates(:,iii),Astates(:,iii),'Color','b');
            else
                plot(Mstates(:,iii),Astates(:,iii),'Color','k');
            end
        end
        
    end
    colormap jet;
    ylabel('Adult Coral');
    xlabel('Macroalgae');
    axis([0 axlim 0 axlim]);
    

    %add nullclines and equilibria
    syms As Ms
    R=-(k(1)*p(1)*(e(1) + As*f(1)*w(1))*(As + Ms - 1))/(a(1) + n(1) + Ms*s(1) + k(1)*p(1)*(e(1) + As*f(1)*w(1)));
    dy1 = a(1)*R + g(1)*As*(1-R-As-Ms) - b(1)*As*Ms - m(1)*As;
    dy2 = s(1)*Ms*(1-R-As-Ms) + s(1)*R*Ms + b(1)*As*Ms - h(1)*Ms - (z(1)*Ms*o(1)*As)/(1+o(1)*As);
    S1=solve(dy1==0, As);
    S2=solve(dy2==0, As);
    plh=fplot(S1(1),[0 1]);set(plh,'color',[0.5 0.5 0.5],'linewidth',2);
    plh=plot(zeros(1,11),0:0.1:1);set(plh,'color','k','linewidth',2);
    plh=fplot(S1(2),[0 1]);set(plh,'color',[0.5 0.5 0.5],'linewidth',2);
    plh=fplot(S2(2),[0 1]);set(plh,'color','k','linewidth',2);
    for jj=1:size(unstabeq,1)
        scatter(unstabeq(jj,1),unstabeq(jj,2),150,'d','k','LineWidth',2.5);
    end
    for jj=1:size(stabeq,1)
        scatter(stabeq(jj,1),stabeq(jj,2),150,'o','filled','k','LineWidth',2.5);%if stable
    end
    ylabel('Adult Coral');
    xlabel('Macroalgae');
    axis([0 axlim 0 axlim]);

    
    %multicolour phase plot
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
    ylabel('Adult Coral');
    xlabel('Macroalgae');
    axis([0 axlim 0 axlim]);
    caxis([0.05 0.1])
    if rep~=length(es)
        clearvars -except es showfigs storestates keepsaddle;
    end
    
end
