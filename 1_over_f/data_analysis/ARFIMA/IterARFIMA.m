function [Betass,BetassErr,Hss,HssErr,Destss,DestErr] = IterARFIMA(N,ds,repeat,stdx)

%{ 
N=100; ds=[0.1:0.1:1]; repeat=100; stdx=20;
%}

% calcualtion
Betass=[];Hss=[];BetassErr=[];HssErr=[];Destss=[];DestErr=[];
for j=1:length(ds)
    d=ds(j);
    Betas=[];Hs=[];
    for i=1:repeat 
        % [Z] = ARFIMA_SIM(N,[],[],d,stdx); % Good H for shorter length
        Z=dgp_arfima(0,[],[],N,stdx,d,0); % Good H for all length
        y=Z+750;
        % PSA
        Fs=1.3; [~,~,beta] = oneoverf(y,Fs);
        % DFA
        [~,H]=DFA_main(y);
        % append the values
        Betas = [Betas beta];
        Hs = [Hs H];  
        % estimate d using Luca's MSE_VARFI
        % (http://www.lucafaes.net/LMSE-MSE_VARFI.html)
        d_min=-0.5; d_max=1;
        Dest=d_estimation(y,d_min,d_max);
    end
Betass=[Betass;Betas];
BetassErr=[BetassErr; std(Betas)/sqrt(N)];
Hss=[Hss;Hs];
HssErr=[HssErr; std(Hs)/sqrt(N)];
Destss=[Destss; Dest];
DestErr=[DestErr; std(Destss)/sqrt(N)];
end
end

