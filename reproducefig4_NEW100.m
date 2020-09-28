scen=1;
relz=50;
ee=[];
thisparams=resultsscen_final(scen).deltacoral(relz).params;
thisparams(:,14)=resultsscen_final(scen).deltacoral(relz).params(:,14)/2;
thisparams=vertcat(thisparams,thisparams);
TP=thisparams;

et=normrnd(0.75,0.01,[100 1]);
%et=thisparams(:,14);
et(et>1)=1; 
et(et<0)=0;
thisparams(:,14)=et;
ee(:,1)=et;
reproducefig3_100 ;

%disrupt=sort(datasample(15:35,5,'Replace',false));
disrupt=5:10:40;
for i=1:length(disrupt)
    thisdl=5;%randi(5)-1;
    et(disrupt(i):(disrupt(i)+thisdl))=normrnd(0.1,0.01,[thisdl+1 1]);
end
et(et>1)=1; 
et(et<0)=0;
ee(:,2)=et;
thisparams(:,14)=et;

reproducefig3_100;  


et=normrnd(0.1,0.01,[100 1]);
%et=thisparams(:,14);
et(et>1)=1; 
et(et<0)=0;
thisparams(:,14)=et;
ee(:,3)=et;
reproducefig3_100 ;
figure;plot(ee(:,1:3))

et=ee(:,1);
disrupt=5:10:80;
for i=1:length(disrupt)
    thisdl=5;%randi(5)-1;
    et(disrupt(i):(disrupt(i)+thisdl))=normrnd(0.1,0.01,[thisdl+1 1]);
end
et(et>1)=1; 
et(et<0)=0;
ee(:,4)=et;
thisparams(:,14)=et;

reproducefig3_100;  

% thisparams=TP;
% % thisparams(:,3)=0;
% % thisparams(:,6)=0;
% % thisparams(:,8)=0;
% % thisparams(:,9)=0;
% thisparams(:,14)=ee(:,1);
% reproducefig3_100_noalgae;  
% 
% thisparams(:,14)=0;
% reproducefig3_100_noalgae; 