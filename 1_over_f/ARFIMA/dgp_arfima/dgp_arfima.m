%--------------------------------------------------------------------------
% Simulation of an ARFIMA(p,d,q) process
% Results=dgp_arfima(c,ARo,MAo,T,SD,d);
%--------------------------------------------------------------------------
% Input:
% c : Constant term
% AR: Parameters of the AR process (column vector: [AR(1),...AR(p)]')
% MA: Parameters of the MA process (column vector: [MA(1),...MA(q)]')
% T : Simulated sample size
% SD: Standard Deviation of the white noise
% d : Fractional memory parameter. d in [0,inf)
% F : Figures; plot of:
%               1.- the series 
%               2.- the sample FAC
%     [possible values: 0-1]
%               1: Figure
% EXAMPLES : 
% ARFIMA(2,0.5,1) >>>>  Results=dgp_arfima(1,[0.5 0.1],[0.2],100,1,0.5);
% ARFIMA(0,0.8,0) >>>>  Results=dgp_arfima(0,[],[],100,1,0.8);
%--------------------------------------------------------------------------
% Output:
% sim           : Simulation of the arfima process 
%--------------------------------------------------------------------------
% Note 1: Follow instructions because this code is not a fool's proof. 
% Note 2: Autocorrelation function is 20% of the sample size. Of course, 
% you can change it.
% Note 2: Fractional memory parameter d cover both stationary and
% nonstationary regions. 
% Note 3: This code uses the Fast Fractional Difference algorithm proposed
% by Jensen & Nielsen (2014)for computational speed purposes.
% Note 4: Below you will find an UAT to check the validity of this code.
% The UAT requires Matlab codes for local Whittle estimation from 
% http://shimotsu.web.fc2.com/Site/Matlab_Codes.html
%--------------------------------------------------------------------------
% Written by: Carlos Vladimir Rodríguez Caballero
% vrodriguez@creates.au.dk
% Aarhus University
% CREATES
%--------------------------------------------------------------------------

function sim=dgp_arfima(c,ARo,MAo,T,SD,d,F);

p=length(ARo); q=length(MAo);
%--------------------------------------------------------------------------
% Simulation of the process
%--------------------------------------------------------------------------
% MA(q) noise:
%--------------------------------------------------------------------------
if SD==0;SD=1;end;
if q==0;
    u=normrnd(0,SD,T+p,1);
else;
    uM=zeros(T,q+1);  
    u=normrnd(0,SD,T+q,1);
    MAb=[1;MAo];
    for j=1:q+1;
        uM(:,j)=MAb(j)*u(q+1-j+1:T+q-j+1,1);
    end;
    u=sum(uM,2);
end;

%--------------------------------------------------------------------------
% d
%--------------------------------------------------------------------------
d=-d;
if d == 0
    fu=u;
else
    fu=fastfrac(u,d);
end

%--------------------------------------------------------------------------
% AR process:
%--------------------------------------------------------------------------
if p==0;
    y=fu;
else;
    yb=zeros(T+p,1);
    for j=p+1:T+p;
        for i=1:p;
            yb(j)=yb(j)+ARo(i)*yb(j-i);
        end;
        yb(j)=c+yb(j)+fu(j-p);
    end;
    y=yb(p+1:T+p,1);
end;
sim=y;

if F==1;
    subplot(2,1,1);
        plot(y);title('Simulated process');
    subplot(2,1,2);
        autocorr(y,length(y)*0.2);  
end;
end

% FAST FRACTIONAL DIFFERENCE ALGORITHM
function[dx] = fastfrac(x,d);
    T = size(x, 1);
    np2 = 2.^nextpow2(2*T - 1);
    k = (1:T-1)';
    b = [1; cumprod((k - d - 1)./k)];
    dx = ifft(fft(x, np2).*fft(b, np2));
    dx = dx(1:T,:);
end

% USER ACCEPTANCE TESTING (UAT)
% for h=1:100
% for d=.1:.1:1.5
% j=round(d*10);
% Z = dgp_arfima(0,[0.1;0.5],[0.5;0.1],5000,1,d,0);
% dmt = dmtrend(Z);
% options = optimset('fminbnd');
% dw(h,j) = fminbnd('whittle',-1,3,options,Z,fix(1000^0.6));
% elw1(h,j) = felw2st(Z,fix(1000^0.6),3,10);					
% elw2(h,j) = felw2st(dmt,fix(1000^0.6),3,10);					
% melw(h,j) = fminbnd('extwhittle',-0.5,2,[],Z,fix(1000^0.6));  
% end
% end
% d=.1:.1:1.5;
% d2=mean(dw);
% d3=mean(elw1);
% d4=mean(elw2);
% d5=mean(melw);
% D=[d',d2',d3',d4',d5']
%     
   
   
