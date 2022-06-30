function [Betass,BetassErr,Hss,HssErr] = IterARFIMA(N,ds,repeat,stdx)

%{ 
N=100; ds=[0.1:0.1:1]; repeat=100; stdx=20;
%}

% calcualtion
Betass=[];Hss=[];BetassErr=[];HssErr=[];
for j=1:length(ds)
    d=ds(j);
    Betas=[];Hs=[];
    for i=1:repeat 
        % [Z] = ARFIMA_SIM(N,[],[],d,stdx); % Good H for shorter length
        Z=dgp_arfima(0,[],[],N,stdx,d,0); % Good for all length
        y=Z+750;
        % PSA
        Fs=1.3; [~,~,beta] = oneoverf(y,Fs);
        % DFA
        [~,H]=DFA_main(y);
        % append the values
        Betas = [Betas beta];
        Hs = [Hs H];  
    end
Betass=[Betass;Betas];
BetassErr=[BetassErr; std(Betas)/sqrt(N)];
Hss=[Hss;Hs];
HssErr=[HssErr; std(Hs)/sqrt(N)];
end
end

