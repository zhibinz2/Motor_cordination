%% Compute errors
% cd /ssd/zhibin/1overf/20220324
% % compute the error
% SynError=[]; % for synchronization and syncopation
% % find the time difference with the closest botton press 
% for i=1:480
%     [minValue,closetIndex]=min(abs(BottonPressTime-PhotocellTime(i))); % closetIndex in BottonPressTime
%     SynError(i)=PhotocellTime(i)-BottonPressTime(closetIndex)
% end
% plot(SynError,'ro');title('synchronization error');xlabel('taps');ylabel('timing error (s)');
% 
% RanError=[]; % for randomization
% % find the time difference with the next botton press
% for i=491:length(PhotocellTime)
%     diffInd=(find((BottonPressTime-PhotocellTime(i))>0));
%     RanError(i)=BottonPressTime(diffInd(1))-PhotocellTime(i);
% end
% RanError=RanError(491:end);
% plot(RanError,'ro'); title('randomization error');xlabel('taps');ylabel('timing error (s)');

cd /ssd/zhibin/1overf/20220331
cd /ssd/zhibin/1overf/20220408
% compute the error for condition 1
clear Error1
Error1=[]; % for synchronization (condition 1)
StartStim1=4+240*(find(allPerm==1)-1)+1; % the starting photocell of condition 1
% find the time difference with the closest botton press 
for i=StartStim1:StartStim1+239 % i=StartStim1
% for i=StartStim1:length(locs) % i=StartStim1
    [minValue,closetIndex]=min(abs(BottonPressTime-PhotocellTime(i))); % closetIndex in BottonPressTime
    Error1(i-StartStim1+1)=BottonPressTime(closetIndex)-PhotocellTime(i);
end
plot(Error1,'r.');title('synchronization error');xlabel('taps');ylabel('timing error (s)');
% remove error > 1 s
indices=find(Error1>1);
Error1(indices)=NaN;
Error1(indices)=[];
% remove error < -1 s
indices=find(Error1<-1);
Error1(indices)=NaN;
Error1(indices)=[];
% Error1(indices)=2;

% compute the error for condition 2
clear Error2
Error2=[]; % for syncopation (condition 2)
StartStim2=4+240*(find(allPerm==2)-1)+1; % the starting photocell of condition 2
% find the time difference with the closest botton press 
for i=StartStim2:StartStim2+238 % i=StartStim2 % syncopation has one press less than condition 1 and 3
    Midpoints=(PhotocellTime(i)+PhotocellTime(i+1))/2;
    [minValue,closetIndex]=min(abs(BottonPressTime-Midpoints)); % closetIndex in BottonPressTime
    Error2(i-StartStim2+1)=BottonPressTime(closetIndex)-Midpoints;
end
plot(Error2,'r.');title('syncopation error');xlabel('taps');ylabel('timing error (s)');
% remove error > 1 s and < -1 s
indices=find(Error2>1);
Error2(indices)=NaN;
Error2(indices)=[];
% Error2(indices)=2;
indices=find(Error2<-1);
Error2(indices)=NaN;
Error2(indices)=[];
% Error2(indices)=2;

% compute the error for condition 3
clear Error3
Error3=[]; % for syncopation (condition 2)
StartStim3=4+240*(find(allPerm==3)-1)+1; % the starting photocell of condition 3
% find the time difference with the next botton press
for i=StartStim3:StartStim3+239 % i=StartStim3 % syncopation has one press less than condition 1 and 3
    diffInd=(find((BottonPressTime-PhotocellTime(i))>0));
    Error3(i-StartStim3+1)=BottonPressTime(diffInd(1))-PhotocellTime(i);
end
plot(Error3,'r.');title('randomization error');xlabel('taps');ylabel('timing error (s)');
% large errors are misses
% remove error > 1 s (In a simple RT task, anything longer than 1 second should be thrown out.)
indices=find(Error3>1);
Error3(indices)=NaN;
Error3(indices)=[];
% Error3(indices)=3;

%% auto correlation and spectrum
% open /ssd/zhibin/1overf/20220324/explore_1_over_f.m

SupTitles={'Synchronization' 'Syncopation' 'RT condition'};

for i=1:3 % i=2;
    if i==1
        y=Error1;
    elseif i==2
        y=Error2;
    else
        y=Error3;
    end
    
    % New FFT
    fcoef=fft(y);
    N=length(y);
    fcoef=fcoef/N;
    halfN=floor(N/2);
    Fs=1/2;% assuming a sampling frequency of 0.5Hz (2 second per sample)
    freqs=linspace(0,Fs/2,halfN+1); % same df
    fcoef=2*fcoef(1:halfN);
    amplitude = abs(fcoef);
    p=amplitude.^2;
    
    % Old FFT 
%     F=fft(y);
%     p1=abs(F).^2; % plot(p1)
%     Fs=1/2;% assuming a sampling frequency of 0.5Hz (2 second per sample)
%     freqs1=linspace(0,Fs/2,length(y)/2+1);

    % plot(log10(freqs),log10(p(1:length(freqs))));
    % xlabel('log of frequency(Hz)');
    % ylabel('log of power');
    % title('{\bf LogLogPlot}');

    % loglog(x,y)
    % grid on


    % 3 plots
    figure('units','normalized','outerposition',[0 0 1 0.3]);
    
    subplot(1,3,1);
    plot(y,'k.');
    xlabel('Trial');ylabel('Error(sec.)');title('A. Time Series');
    
    subplot(1,3,2);
    autocorr(y,'NumLags',200,'NumSTD',2);
    % [acf,lags,bounds] = autocorr(RTs(ConditionExplore),'NumLags',25,'NumSTD',2);
    xlabel('Lag');ylabel('Corr.');title('B. Autocorr.');
    ylim([-0.4 1]);
    
    subplot(1,3,3);
    %xx=log10(freqs(2:end));yy=log10(p(2:length(freqs)));
    xx=log10(freqs(2:end));yy=log10(p(1:length(p)));
    plot(xx,yy);hold on;
    tbl=table(xx',yy');
    mdl=fitlm(tbl,'linear');
    plotAdded(mdl);
    xlabel('Log(f)');ylabel('Log(power)');title('C. Spectrum');
    %ylim([-3 3]);
    ylim([-7 0]);
    
    suptitle([char(SupTitles(i)) ' in subject ' num2str(seed)]);
    
%     % save figure handles
%     ax(i,1)=gca;
 
%     % save as .fig file
%     saveas(gcf,[char(SupTitles(i)) '.jpg'])
end

% % Put all 3 picutre together
% h4 = figure; %create new figure
% s1 = subplot(3,1,1); %create and get handle to the subplot axes
% s2 = subplot(3,1,2);
% s3 = subplot(3,1,3);
% fig1 = get(ax(1),'children'); %get handle to all the children in the figure
% fig2 = get(ax(2),'children');
% fig3 = get(ax(3),'children');
% copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
% copyobj(fig2,s2);
% copyobj(fig3,s3);

% close all;
% fig(1) = openfig('Synchronization.fig'); % open figure
% fig(2) = openfig('Syncopation.fig'); % open figure
% fig(3) = openfig('Randomization.fig'); % open figure
% new_fig = figure;
% ax_new = gobjects(size(fig));
% for i=1:3
%     ax = subplot(3,1,i);
%     ax_old = findobj(fig(i), 'type', 'axes');
%     ax_new(i) = copyobj(ax_old,new_fig);
%     ax_new(i).YLimMode = 'manual';
%     ax_new(i).Position = ax.Position;
% %     ax_new(i).Position(4) = ax_new(i).Position(4)-0.02;
%     delete(ax);
% end