Error1
Error2
Error3

%% matlab example 1
% https://www.mathworks.com/matlabcentral/fileexchange/19795-detrended-fluctuation-analysis

% DATA=Error1;
% win_length=100;
% order=1;

% function output1=DFA(DATA,win_length,order)
%       
% N=length(DATA);   % length of the data
% n=floor(N/win_length); % the number of windows
% N1=n*win_length; % length ot the data again
% 
% y=zeros(N1,1); % initialize y 
% Yn=zeros(N1,1);
% 
% fitcoef=zeros(n,order+1); % initialize fitcoef
% 
% mean1=mean(DATA(1:N1)); % grand mean of the data
% 
% % sum of deviations from the grand mean in all data length
% for i=1:N1      
%     y(i)=sum(DATA(1:i)-mean1); 
% end
% y=y';
% 
% % % my alternative method, use cumsum instead (but incorrect, this will result in original length)
% % y=cumsum(DATA-mean1);
% 
% for j=1:n
%     fitcoef(j,:)=polyfit(1:win_length,y(((j-1)*win_length+1):j*win_length),order);
% end
% %  [P] = POLYFIT(X,Y,N) returns the polynomial coefficients P for use with POLYVAL to obtain error estimates for
% %   predictions.  
% 
% for j=1:n
%     Yn(((j-1)*win_length+1):j*win_length)=polyval(fitcoef(j,:),1:win_length);
% end
% %   Y = POLYVAL(P,X) returns the value of a polynomial P evaluated at X. P
% %   is a vector of length N+1 whose elements are the coefficients of the
% %   polynomial in descending powers: Y = P(1)*X^N + P(2)*X^(N-1) + ... + P(N)*X + P(N+1)
% 
% sum1=sum((y'-Yn).^2)/N1;
% sum1=sqrt(sum1);
% output1=sum1;

%% https://www.mathworks.com/matlabcentral/fileexchange/19795-detrended-fluctuation-analysis
DATA=repmat(Error1,1,5);
DATA=DATA(1:512); % use only the first 512 minimal sample number, according to Marmelat & Deligniere 2012

DATA=Error1;
DATA=Error2;
DATA=Error3;

function [D,Alpha1]=DFA_main(DATA)
% DATA should be a time series of length(DATA) greater than 2000,and of column vector.
%A is the alpha in the paper
%D is the dimension of the time series % ?
%n can be changed to your interest

% n=100:100:1000;
% n=2:2:200;
n=2:2:512
N1=length(n);
F_n=zeros(N1,1);
 for i=1:N1
     F_n(i)=DFA(DATA,n(i),1);
 end
 n=n';
 plot(n,F_n);
%  plot(log(n),log(F_n));
xlabel('n')
ylabel('F(n)')
 A=polyfit(log(n(1:end)),log(F_n(1:end)),1);
Alpha1=A(1);
 D=3-A(1);%  ?
return



%% https://www.mathworks.com/matlabcentral/fileexchange/60026-r-dfa-robust-detrended-fluctuation-analysis
% functions in folder "r-DFA"

%% Mike Cohen lecture (https://www.youtube.com/watch?v=-RmxLZF8adI)
