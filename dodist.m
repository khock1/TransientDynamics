% Script to implement distrubances for the coral  model

% This code supplements the article "Transient dynamics mask the resilience
% of coral reefs" by Hock et al.
% Author & copyright: Karlo Hock, University of Queensland. 2021

%only some functionalities are used in the paper

function [pt,at,st,nt,gt,bt,mt,ht,zt,ot,wt,ft,kt,et]=dodist(p,a,s,n,g,b,m,h,z,o,w,f,k,e,timesteps,cyctimes,blchtimes,graztimes,conntimes, grazrecovery, morttimes, hmorttimes, cgrowvartimes, spcAgrowtimes)


at = repmat(a,[timesteps 1]);%0.35
st = repmat(s,[timesteps 1]);%0.8
nt = repmat(n,[timesteps 1]);%0.45
gt = repmat(g,[timesteps 1]);%0.2 0.4
bt = repmat(b,[timesteps 1]);%0.1
mt = repmat(m,[timesteps 1]);%0.05 078
ht = repmat(h,[timesteps 1]);%0.2 0.01
zt = repmat(z,[timesteps 1]);%0.64 0.28
ot = repmat(o,[timesteps 1]);%4
wt = repmat(w,[timesteps 1]);%0.4
ft = repmat(f,[timesteps 1]);%1
pt = repmat(p,[timesteps 1]);%1
kt = repmat(k,[timesteps 1]);%0.1
et = repmat(e,[timesteps 1]);


for i=1:timesteps
    wn(i)=-0.204*i+1;
end

for t=1:timesteps
    if ismember(t,morttimes)
        mt(t)=awgn(m,5,'measured');
    end

    if ismember(t,graztimes)
        zt(t)=0.64*rand;%0.05
    else
        zt(t)=z;
    end
    if ismember(t,conntimes)
        et(t)=awgn(e,5,'measured');%this

    else

        et(t)=et(t);
    end
    if ismember(t,hmorttimes)
        ht(t)=awgn(h,5,'measured');
    end

    if ismember(t,cgrowvartimes)
        gt(t)=awgn(g,5,'measured');
    end

    if ismember(t,spcAgrowtimes)
        st(t)=awgn(s,5,'measured');
    end
    if ismember(t,cyctimes)

        mt(t)=0.5*rand;
        ht(t)=0.8*rand;
    else
        mt(t)=mt(t);
        ht(t)=ht(t);
    end
    if ismember(t,blchtimes)
        mt(t)=0.8*rand;
    else
        mt(t)=mt(t);
    end
    mt(mt>1)=1;
    mt(mt<0)=0;
    st(st>1)=1;
    st(st<0)=0;
    ht(ht>1)=1;
    ht(ht<0)=0;
    gt(gt>1)=1;
    gt(gt<0)=0;
    et(et>1)=1;
    et(et<0)=0;

end

if grazrecovery < timesteps
    incr=(z-0.05)/6;
    incr0=0.05+incr;
    for t=1:6
        zt(grazrecovery+1+(t-1))=incr0;
        incr0=incr0+incr;
    end
end

end

