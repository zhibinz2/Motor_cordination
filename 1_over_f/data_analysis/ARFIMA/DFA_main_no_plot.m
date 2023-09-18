function [D,Alpha1,n,F_n,FitValues]=DFA_main_no_plot(DATA)
% DATA should be a time series of length(DATA) greater than 2000,and of column vector.
% But aslo works for 100 samples in my data.
% A is the alpha in the paper
% (https://aip.scitation.org/doi/pdf/10.1063/1.166141)
% D is the dimension of the time series % ?
% n can be changed to your interest

% n=100:100:1000; % win_lengths (function default)
% n=2:2:200; % all win_lengths(customized to our need)
n=5:5:length(DATA)/2; % all 32 different win_lengths (customized to our need)
Nw=length(n); % number of different win_lengths
F_n=zeros(Nw,1); % initialize RMS values in the FDA time series
 for i=1:Nw
     F_n(i)=DFA(DATA,n(i),1);
 end
n=n';
 

% A=polyfit(log10(n(scaleInds)),log10(F_n(scaleInds)),1);
A=polyfit(log10(n),log10(F_n),1);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
D=3-A(1); % ? D is the dimension of the time series % ?

return
