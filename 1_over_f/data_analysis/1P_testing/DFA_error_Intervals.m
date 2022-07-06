%% plot DFA on errors
clear
cd /ssd/zhibin/1overf/20220331/
cd cd /ssd/zhibin/1overf/20220408
load /ssd/zhibin/1overf/20220331/clean_0425/20220331workspace.mat
Error1; Error2; Error3;
figure;
for condi=1:3
    if condi==1;y=Error1;end
    if condi==2;y=Error2;end
    if condi==3;y=Error3;end
    subplot(1,3,condi);
    [D,Alpha1,n,F_n,FitValues]=DFA_main(y);
%     % figure; plot(Error1);
%     figure;
    plot(log10(n),log10(F_n),'bx');
    xlabel('Scale [log10(n)]') % win_lengths 
    ylabel('RMS [log10(F-n)]') % RMS values
    title(['DFA condi ' num2str(condi)]);
    hold on;
    plot(log10(n),FitValues,'r--');
    legend({'Data',['Fit (DFA=' num2str(Alpha1) ')']},'Location','southeast');
    ylim([-1.4 0.4]);
end
suptitle(['DFA of errors' num2str(seed)]);

%% plot DFA on intervals

% compute intervals
Intv1=[];Intv2=[];Intv3=[];
% condition 1
StartStim1=4+240*(find(allPerm==1)-1)+1; % the starting photocell of condition 1
% find the time difference with the closest botton press 
[minValue1,closetIndex1]=min(abs(BottonPressTime-PhotocellTime(StartStim1))); % closetIndex in BottonPressTime
[minValue2,closetIndex2]=min(abs(BottonPressTime-PhotocellTime(StartStim1+239))); 
% [minValue2,closetIndex2]=min(abs(BottonPressTime-PhotocellTime(end))); 
BP01=BottonPresTime01(BottonPresTimeInd(closetIndex1):BottonPresTimeInd(closetIndex2)); % plot(BP01)
Intv1=Calinterval(BP01')/2000;% plot(BPint) % in seconds
% condition 2
StartStim2=4+240*(find(allPerm==2)-1)+1; % the starting photocell of condition 2
% find the time difference with the closest botton press 
[minValue1,closetIndex1]=min(abs(BottonPressTime-PhotocellTime(StartStim2))); % closetIndex in BottonPressTime
[minValue2,closetIndex2]=min(abs(BottonPressTime-PhotocellTime(StartStim2+239))); 
BP01=BottonPresTime01(BottonPresTimeInd(closetIndex1):BottonPresTimeInd(closetIndex2)); % plot(BP01)
Intv2=Calinterval(BP01')/2000;% plot(Intv2) % in seconds
% condition 3
StartStim3=4+240*(find(allPerm==3)-1)+1; % the starting photocell of condition 3
% find the time difference with the next botton press
[minValue1,closetIndex1]=min(abs(BottonPressTime-PhotocellTime(StartStim3))); % closetIndex in BottonPressTime
[minValue2,closetIndex2]=min(abs(BottonPressTime-PhotocellTime(StartStim3+239))); 
[minValue2,closetIndex2]=min(abs(BottonPressTime-PhotocellTime(end))); 
BP01=BottonPresTime01(BottonPresTimeInd(closetIndex1):BottonPresTimeInd(closetIndex2)); % plot(BP01)
Intv3=Calinterval(BP01')/2000;% plot(Intv3) % in seconds

% plot DFA
figure;
for condi=1:3
    if condi==1;y=Intv1;end
    if condi==2;y=Intv2;end
    if condi==3;y=Intv3;end
    subplot(1,3,condi);
    [D,Alpha1,n,F_n,FitValues]=DFA_main(y);
    %     % figure; plot(Error1);
    %     figure;
    plot(log10(n),log10(F_n),'bx');
    xlabel('Scale [log10(n)]') % win_lengths 
    ylabel('RMS [log10(F-n)]') % RMS values
    title(['DFA condi ' num2str(condi)]);
    hold on;
    plot(log10(n),FitValues,'r--');
    legend({'Data',['Fit (DFA=' num2str(Alpha1) ')']},'Location','southeast');
    ylim([-1.4 0.4]);
end

suptitle(['DFA of tapping intervals' num2str(seed)]);
