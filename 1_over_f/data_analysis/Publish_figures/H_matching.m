clear

%% addpath
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA

% % d estimate and removal
% % http://www.lucafaes.net/LMSE-MSE_VARFI.html
% addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/MSE-VARFI
% % Granger Causality
% addpath(genpath('/home/zhibinz2/Documents/GitHub/MVGC1'));
% cd /home/zhibinz2/Documents/GitHub/MVGC1
% run startup
% addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/VAR_Granger


%% Color Scheme
% Plots - color scheme
red   = [1 0 0];
pink  = [1 0.65 0.75];
black = [0 0 0];
white = [1 1 1];
blue  = [0 0 1];
mediumblue = [0 0.4 0.7];
green = [0 1 0];
darkgreen = [0 0.5 0];
grey  = [0.5 0.5 0.5];
yellow  = [1 1 0];
deepyellow  = [1 0.8 0.2];
gold = [212/255 175/255 55/255];
brown = [150/255 75/255 0];
megenta = [1 0 1];% fill([0 1 1 0],[0 0 1 1],megenta)
cyan = [0 1 1]; % fill([0 1 1 0],[0 0 1 1],cc)
purple = [0.6 0.1 0.9];
% https://www.mathworks.com/help/matlab/creating_plots/specify-plot-colors.html
matlab_blue=[0 0.4470 0.7410];
matlab_orange=[0.8500 0.3250 0.0980];
matlab_gold=[0.9290 0.6940 0.1250];
matlab_purple=[0.4940 0.1840 0.5560];
matlab_green=[0.4660 0.6740 0.1880];
matlab_cyan=[0.3010 0.7450 0.9330];
matlab_red=[0.6350 0.0780 0.1840];
% combine colors
condicolors=[darkgreen;red;blue;megenta;purple;purple];
dire3colors=[darkgreen;brown;megenta];
% syn2colors=[matlab_blue;matlab_orange];
syn2colors=[darkgreen;pink];
HNLcolors = [darkgreen; deepyellow; pink];

% % test color
% showcolor=black;
% imagesc(cat(3,showcolor(1),showcolor(2),showcolor(3)));

% colormap
cd /home/zhibinz2/Documents/GitHub/matlab-archive/hnlcode/common/gen_code/color
addpath /home/zhibinz2/Documents/GitHub/matlab-archive/hnlcode/common/gen_code/color
hnc = hotncold(100);
hottermap=hotter(100);
% 
% % test colormap
% mycolormap=hottermap;
% mymat=repmat(-0.99:0.02:0.99,100,1); % transitional color
% % mymat=(rand(100)-0.5)*2; % random from -1 to 1
% imagesc(mymat);colorbar;
% colormap(mycolormap);

%% H compute and extract
% SECT 10 redo H and intervals_H_removed and save them separately (matched int)
% in original order
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

cd /ssd/zhibin/1overf/
H_all=[];
for s=1:numSes
    clear intervals
    runid=num2str(seeds(s,:));
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat'],'intervals')   
    for b=1:12
        [~,H_all(1,s,b)]=DFA_main(intervals{b}(:,1));
        [~,H_all(2,s,b)]=DFA_main(intervals{b}(:,2));
    end
end
%%
% H-int matching for L and R as in Marmelat 2012
% organize H_all 
H_all; % 2x12numSesx12trials
H_Lall=squeeze(H_all(1,:,:)); % 12session by 12trials
H_Lall=reshape(H_Lall',[],1); % in time sequence
H_Rall=squeeze(H_all(2,:,:)); % 12session by 12trials
H_Rall=reshape(H_Rall',[],1); % in time sequence
H_all_LR=[H_Lall; H_Rall];


%%
condition_all=[];
for s=1:numSes
    clear conditions
    runid=num2str(seeds(s,:));
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat'],'conditions');
    condition_all(s,:)=conditions; % 12 session x 12 trials
end
% reshape into a vector in time sequence
condition_all=reshape(condition_all',[],1); % 144 x 1 


%%
% indices for 4 states
states4names={'Independent','Leading','Following','Bidrectional'};
% find the indices for each condition in L R conbined sequence (4 states)
% and organize for PLS
uncoupleInd_LR=[find(condition_all==1);12*numSes+find(condition_all==1)];
leadingInd_LR=[find(condition_all==2);12*numSes+find(condition_all==3)];
followingInd_LR=[find(condition_all==3);12*numSes+find(condition_all==2)];
mutualInd_LR=[find(condition_all==4);12*numSes+find(condition_all==4)];
Inds4_LR=[uncoupleInd_LR leadingInd_LR followingInd_LR mutualInd_LR];


%%
% H-int matching for L and R as in Marmelat 2012
% indicies in the synch and synco time sequence
synchind=[1:3 7:9 13:15 19:21 25:27 31:33]; % 3 trials x 6 sessions
syncoind=[4:6 10:12 16:18 22:24 28:30 34:36]; % 3 trials x 6 sessions
synind=[synchind;syncoind];

%% Bonanza plot H matching
% PLOT 10-2 H matching as in Marmelat 2012 (matched int) Reproduced the plot in poster
for subplots3 =1 
% figure; 
canvas(0.5,0.35);
subplot(1,3,1); %##############
for subplot1 =1;
plot(H_Lall(uncoupleInd_LR(synchind(1:3))),H_Rall(uncoupleInd_LR(synchind(1:3))),'.','MarkerSize',30,'color',darkgreen);hold on;
plot(H_Lall(uncoupleInd_LR(synchind(4:6))),H_Rall(uncoupleInd_LR(synchind(4:6))),'square','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(uncoupleInd_LR(synchind(7:9))),H_Rall(uncoupleInd_LR(synchind(7:9))),'^','MarkerSize',8,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(uncoupleInd_LR(synchind(10:12))),H_Rall(uncoupleInd_LR(synchind(10:12))),'pentagram','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(uncoupleInd_LR(synchind(13:15))),H_Rall(uncoupleInd_LR(synchind(13:15))),'*','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(uncoupleInd_LR(synchind(16:18))),H_Rall(uncoupleInd_LR(synchind(16:18))),'diamond','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(uncoupleInd_LR(syncoind(1:3))),H_Rall(uncoupleInd_LR(syncoind(1:3))),'.','MarkerSize',30,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(uncoupleInd_LR(syncoind(4:6))),H_Rall(uncoupleInd_LR(syncoind(4:6))),'square','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(uncoupleInd_LR(syncoind(7:9))),H_Rall(uncoupleInd_LR(syncoind(7:9))),'^','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(uncoupleInd_LR(syncoind(10:12))),H_Rall(uncoupleInd_LR(syncoind(10:12))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(uncoupleInd_LR(syncoind(13:15))),H_Rall(uncoupleInd_LR(syncoind(13:15))),'*','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(uncoupleInd_LR(syncoind(16:18))),H_Rall(uncoupleInd_LR(syncoind(16:18))),'diamond','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
ax = gca;
ax.FontSize = 13; 
xlabel('DFA Hurst Exponent, Participant A','FontSize',15);
ylabel('DFA Husrt Exponent, Participant B','FontSize',15);
title('Independent','FontSize',15);


% A=[];S=[]; Alpha1=[];FitValues=[];
% A=polyfit(H_Lall(uncoupleInd_LR(1:36)),H_Rall(uncoupleInd_LR_LR(1:36)),1);
% Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
% FitValues=polyval(A,H_Lall(uncoupleInd_LR(1:36)));
% hold on; plot(H_Lall(uncoupleInd_LR(1:36)),FitValues,'k-','LineWidth',3);

% fit the synch data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(uncoupleInd_LR(1:2:35)),H_Rall(uncoupleInd_LR(1:2:35)),0.2,1.4);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'-','color',darkgreen,'LineWidth',3);
text(max(xxx),min(FitValues),sprintf('\\rho=%.2f',RHO),'Color',darkgreen,'FontSize', 15)

% fit the synco data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(uncoupleInd_LR(2:2:36)),H_Rall(uncoupleInd_LR(2:2:36)),0.2,1.4);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'-','color',pink,'LineWidth',3);
% xlim([0.5 1]);ylim([0.5 1]);
xlim([0.2 1.4]);ylim([0.2 1.4]);
plot([0.5 1], [0.5 1],'m--');
hold off;
text(max(xxx),max(FitValues),sprintf('\\rho=%.2f',RHO),'Color',pink,'FontSize', 15)
% legend('synch','','','','synco','','','','location','northwest');
% subtitle(['(uncouple: linear regression slope = ' num2str(Alpha1,'%.3f') ')'])
grid on;
end
subplot(1,3,2); %##############
for subplot2 =1;
plot(H_Lall(leadingInd_LR(synchind(1:3))),H_Rall(leadingInd_LR(synchind(1:3))),'.','MarkerSize',30,'color',darkgreen); hold on;
plot(H_Lall(leadingInd_LR(synchind(4:6))),H_Rall(leadingInd_LR(synchind(4:6))),'square','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(leadingInd_LR(synchind(7:9))),H_Rall(leadingInd_LR(synchind(7:9))),'^','MarkerSize',8,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(leadingInd_LR(synchind(10:12))),H_Rall(leadingInd_LR(synchind(10:12))),'pentagram','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(leadingInd_LR(synchind(13:15))),H_Rall(leadingInd_LR(synchind(13:15))),'*','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(leadingInd_LR(synchind(16:18))),H_Rall(leadingInd_LR(synchind(16:18))),'diamond','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(leadingInd_LR(syncoind(1:3))),H_Rall(leadingInd_LR(syncoind(1:3))),'.','MarkerSize',30,'color',pink);
plot(H_Lall(leadingInd_LR(syncoind(4:6))),H_Rall(leadingInd_LR(syncoind(4:6))),'square','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(leadingInd_LR(syncoind(7:9))),H_Rall(leadingInd_LR(syncoind(7:9))),'^','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(leadingInd_LR(syncoind(10:12))),H_Rall(leadingInd_LR(syncoind(10:12))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(leadingInd_LR(syncoind(13:15))),H_Rall(leadingInd_LR(syncoind(13:15))),'*','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(leadingInd_LR(syncoind(16:18))),H_Rall(leadingInd_LR(syncoind(16:18))),'diamond','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Rall(followingInd_LR(synchind(1:3))),H_Lall(followingInd_LR(synchind(1:3))),'.','MarkerSize',30,'color',darkgreen);
plot(H_Rall(followingInd_LR(synchind(4:6))),H_Lall(followingInd_LR(synchind(4:6))),'square','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Rall(followingInd_LR(synchind(7:9))),H_Lall(followingInd_LR(synchind(7:9))),'^','MarkerSize',8,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Rall(followingInd_LR(synchind(10:12))),H_Lall(followingInd_LR(synchind(10:12))),'pentagram','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Rall(followingInd_LR(synchind(13:15))),H_Lall(followingInd_LR(synchind(13:15))),'*','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Rall(followingInd_LR(synchind(16:18))),H_Lall(followingInd_LR(synchind(16:18))),'diamond','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Rall(followingInd_LR(syncoind(1:3))),H_Lall(followingInd_LR(syncoind(1:3))),'.','MarkerSize',30,'color',pink);
plot(H_Rall(followingInd_LR(syncoind(4:6))),H_Lall(followingInd_LR(syncoind(4:6))),'square','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Rall(followingInd_LR(syncoind(7:9))),H_Lall(followingInd_LR(syncoind(7:9))),'^','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(H_Rall(followingInd_LR(syncoind(10:12))),H_Lall(followingInd_LR(syncoind(10:12))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Rall(followingInd_LR(syncoind(13:15))),H_Lall(followingInd_LR(syncoind(13:15))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Rall(followingInd_LR(syncoind(16:18))),H_Lall(followingInd_LR(syncoind(16:18))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
ax = gca;
ax.FontSize = 13; 
xlabel('DFA Hurst Exponent, Leader','FontSize',15);
ylabel('DFA Hurst Exponent, Follower','FontSize',15);
title('Unidirectional','FontSize',15);

% A=[];Alpha1=[];FitValues=[];
% A=polyfit([H_Lall(leadingInd_LR(1:36)); H_Rall(followingInd_LR(1:36))], [H_Rall(leadingInd_LR(1:36)); H_Lall(followingInd_LR(1:36))],1);
% Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
% FitValues=polyval(A,[H_Lall(leadingInd_LR(1:36)); H_Rall(followingInd_LR(1:36))]);
% hold on; plot([H_Lall(leadingInd_LR(1:36)); H_Rall(followingInd_LR(1:36))],FitValues,'k-','LineWidth',3);

% fit the synch data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues([H_Lall(leadingInd_LR(1:2:35)); H_Rall(followingInd_LR(1:2:35))], ...
    [H_Rall(leadingInd_LR(1:2:35)); H_Lall(followingInd_LR(1:2:35))],0.2,1.4);
A=polyfit(xxx, yyy,1);
[RHO,~]=corr(xxx, yyy);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,xxx);
hold on; 
plot(xxx,FitValues,'-','color',darkgreen,'LineWidth',3);
text(max(xxx)-0.1,max(FitValues)-0.03,sprintf('\\rho=%.2f',RHO),'Color',darkgreen,'FontSize', 15)

% fit the synco data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues([H_Lall(leadingInd_LR(2:2:36)); H_Rall(followingInd_LR(2:2:36))],...
    [H_Rall(leadingInd_LR(2:2:36)); H_Lall(followingInd_LR(2:2:36))],0.2,1.4);
A=polyfit(xxx, yyy,1);
[RHO,~]=corr(xxx, yyy);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'-','color',pink,'LineWidth',3);
text(max(xxx)-0.1,max(FitValues)+0.1,sprintf('\\rho=%.2f',RHO),'Color',pink,'FontSize', 15)
% xlim([0.5 1]);ylim([0.5 1]);
xlim([0.2 1.4]);ylim([0.2 1.4]);
plot([0.5 1], [0.5 1],'m--');
plot([0.5 1], [0.5 1],'m--');
hold off;
% legend('synch','','','','synco','','','','','','','','','','','','location','northwest');
% subtitle(['(leading or following: linear regression slope = ' num2str(Alpha1,'%.3f') ')'])
grid on;
end
subplot(1,3,3); %##############
for subplot3 =1;
plot(H_Lall(mutualInd_LR(synchind(1:3))),H_Rall(mutualInd_LR(synchind(1:3))),'.','MarkerSize',30,'color',darkgreen);hold on;
plot(H_Lall(mutualInd_LR(synchind(4:6))),H_Rall(mutualInd_LR(synchind(4:6))),'square','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(mutualInd_LR(synchind(7:9))),H_Rall(mutualInd_LR(synchind(7:9))),'^','MarkerSize',8,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(mutualInd_LR(synchind(10:12))),H_Rall(mutualInd_LR(synchind(10:12))),'pentagram','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(mutualInd_LR(synchind(13:15))),H_Rall(mutualInd_LR(synchind(13:15))),'*','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(mutualInd_LR(synchind(16:18))),H_Rall(mutualInd_LR(synchind(16:18))),'diamond','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(mutualInd_LR(syncoind(1:3))),H_Rall(mutualInd_LR(syncoind(1:3))),'.','MarkerSize',30,'color',pink);
plot(H_Lall(mutualInd_LR(syncoind(4:6))),H_Rall(mutualInd_LR(syncoind(4:6))),'square','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(mutualInd_LR(syncoind(7:9))),H_Rall(mutualInd_LR(syncoind(7:9))),'^','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(mutualInd_LR(syncoind(10:12))),H_Rall(mutualInd_LR(syncoind(10:12))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(mutualInd_LR(syncoind(13:15))),H_Rall(mutualInd_LR(syncoind(13:15))),'*','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(mutualInd_LR(syncoind(16:18))),H_Rall(mutualInd_LR(syncoind(16:18))),'diamond','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
ax = gca;
ax.FontSize = 13; 
xlabel('DFA Hurst Exponent, Participant A','FontSize',15);
ylabel('DFA Husrt Exponent, Participant B','FontSize',15);
title('Bidirectional','FontSize',15);
% xlim([0.5 1]);ylim([0.5 1]);
xlim([0.2 1.4]);ylim([0.2 1.4]);
% fit all the data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(mutualInd_LR(1:36)),...
    H_Rall(mutualInd_LR(1:36)),0.2,1.4);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1); 
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'k-','LineWidth',3);
text(max(xxx)-0.2,max(FitValues)-0.2,sprintf('\\rho=%.2f',RHO),'Color',[0 0 0],'FontSize', 15)

% fit the synch data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(mutualInd_LR(1:2:35)),...
    H_Rall(mutualInd_LR(1:2:35)),0.2,1.4);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'-','color',darkgreen,'LineWidth',3);
text(max(xxx)-0.1,max(FitValues)-0.1,sprintf('\\rho=%.2f',RHO),'Color',darkgreen,'FontSize', 15)

% fit the synco data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(mutualInd_LR(2:2:36)),...
    H_Rall(mutualInd_LR(2:2:36)),0.2,1.4);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'-','color',pink,'LineWidth',3);
text(max(xxx)-0.5,max(FitValues)-0.25,sprintf('\\rho=%.2f',RHO),'Color',pink,'FontSize', 15)
% xlim([0.5 1]);ylim([0.5 1]);
xlim([0.2 1.4]);ylim([0.2 1.4]);
plot([0.5 1], [0.5 1],'m--');
hold off;
% subtitle(['(mutual): linear regression slope = ' num2str(Alpha1,'%.3f') ')'])
grid on;
end
lg=legend('synch - subj pair 1','synch - subj pair 2','synch - subj pair 3','synch - subj pair 4', ...
    'synch - subj pair 5','synch - subj pair 6', ...
    'synco - subj pair 1','synco - subj pair 2','synco - subj pair 3','synco - subj pair 4',...
    'synco - subj pair 5','synco - subj pair 6','location','eastoutside');
lg.Position = [0.9475 0.15 0.01 0.75];
% sgtitle(['H matching as in Marmelat 2012  ^{* PLOT 10-2}']);
delete(findall(gcf,'type','annotation'))
sg=annotation('textbox',[0.07 0.01 0.22 0.05],'string',...
    'H matching as in Marmelat 2012  ^{* PLOT 10-2}')
sg.Rotation=90
set(gcf,'color','w'); % set background white for copying in ubuntu
end



