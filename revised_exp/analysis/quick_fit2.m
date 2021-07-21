%quici.fit for kenta
%simple line fit for d'
clear all 
close all
clc

addpath ..\data;
subcode = 'kw2';

nstart=2;      %start on 6th trial
dprange=[0.0001 5];
targetdprime=2;
S = load([subcode 'data.mat']);

nblk = length(S.data);

siglev = [];
iresp=[];
icor=[];

for kk = 1:nblk
    siglev = [siglev; S.data{kk}.track(nstart:end,1)];
    icor = [icor; S.data{kk}.track(nstart:end,2)];
    iresp = [iresp; S.data{kk}.track(nstart:end,3)];
end
sigleva=round(siglev*2000);
sigleva=sigleva/2000;
correctness=zeros(size(iresp));
indices=find(iresp==icor);
correctness(indices)=1;

levvals=unique(sigleva);
nlevs=length(levvals);
%this is number of levels tested. 
newlev=[]; nhit=[]; nfa=[]; nsig=[]; nnosig=[];
dprime=[];
for kk=1:nlevs
    newlev(kk)=levvals(kk);
    idx=find(sigleva==newlev(kk));
    ticor=icor(idx);
    tiresp=iresp(idx);
    pc(kk)=length(find(tiresp==ticor))/length(ticor);
    tidx=find(ticor==1 & tiresp==1);
    nhit(kk)=length(tidx);
    nsig(kk)=length(find(ticor==1));
    hit(kk)=nhit(kk)/nsig(kk);
    tidx=find(ticor==2 & tiresp==1);
    nfa(kk)=length(tidx);
    nnosig(kk)=length(find(ticor==2));
    fa(kk)=nfa(kk)/nnosig(kk);
    dprime(kk)=norminv(hit(kk))-norminv(fa(kk));
    ntrials(kk)=length(idx);
end

idx=find((dprime>=dprange(1)) & (dprime<=dprange(2)));
levvals=levvals(idx);
dprime=dprime(idx);
pc=pc(idx);
ntrials=ntrials(idx);

[b,se_b,mse]=lscov([ones(size(levvals)) levvals],dprime',ntrials);
tottrls=sum(ntrials);
thres_est=(targetdprime-b(1))/b(2);
predline=levvals*b(2)+b(1);

% [b2,se_b2,mse2]=lscov([ones(size(levvals)) levvals],dprime');
% pred2=levvals*b2(2)+b2(1);  xx this is without weighting, which I removed
% plot(levvals,dprime,'k-',levvals,predline,'k--',thres_est,targetdprime,'or',...
%     levvals,pred2,'b--');

plot(levvals,dprime,'k-',levvals,predline,'k--',thres_est,targetdprime,'or');
hold
scatter(levvals,dprime,5*ntrials,'k');
hold off
xlabel('level as % deltaf/f');
ylabel('dprime');
title(['obs ' subcode '   total trials ' int2str(tottrls) '    thresh est '...
    num2str(thres_est)]);
legend('dprime','fit','thres est');
% legend('dprime','fit','thres est','unweighted');

[levvals ntrials']

%here a problem is there are so few values that most get tossed.
%this needs to be addressed.
    