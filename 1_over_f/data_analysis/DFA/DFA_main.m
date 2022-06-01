function [D,Alpha1]=DFA_main(DATA)
% DATA should be a time series of length(DATA) greater than 2000,and of column vector.
% A is the alpha in the paper
% (https://aip.scitation.org/doi/pdf/10.1063/1.166141)
% D is the dimension of the time series % ?
% n can be changed to your interest

% n=100:100:1000; % win_lengths (function default)
% n=2:2:200; % all win_lengths(customized to our need)
n=10:30:300; % all 32 different win_lengths (customized to our need)
Nw=length(n); % number of different win_lengths
F_n=zeros(Nw,1); % initialize RMS values in the FDA time series
 for i=1:Nw
     F_n(i)=DFA(DATA,n(i),1);
 end
 n=n';
 
%  plot(n,F_n);
plot(log(n),log(F_n));
xlabel('Scale - log(n)') % win_lengths 
ylabel('RMS - log(F-n)') % RMS values
title('DFA')

A=polyfit(log(n(1:end)),log(F_n(1:end)),1);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
D=3-A(1); % ? D is the dimension of the time series % ?

% plot the fit
FitValues=polyval(A,log(n(1:end)));
hold on;
plot(log(n(1:end)),FitValues,'r--');
legend({'Data',['Fit (DFA=' num2str(Alpha1) ')']});
return
