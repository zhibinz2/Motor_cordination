stimTrials_synchronization_r,
stimTrials_syncopation_r,
stimTrials_rt_r,
tapTrials_synchronization_r,
tapTrials_syncopation_r,
tapTrials_rt_r

%% select which two data set to process
stimTrial_data=stimTrials_synchronization_r;
tapTrial_data=tapTrials_synchronization_r;

%% low pass filter < 8 Hz
% add paddings
padding=zeros(round(size(stimTrial_data,1)/2), size(stimTrial_data,2),size(stimTrial_data,3));
stimTrial_paddings=cat(1,padding,stimTrial_data,padding);
close;figure;
for i=1:16
    plot(stimTrial_paddings(:,:,i));pause(0.5);
end

% which -all filtfilt.m
rmpath /home/zhibin/Documents/GitHub/matlab-archive/hnlcode/external/fieldtrip-20141007
cd  /home/zhibin/Documents/GitHub/matlab-archive
cd /usr/local/MATLAB/R2019a/toolbox/signal/signal
which filtfilt.m

% high pass (no paddings needed)
Hd = makefilter(sr,1,0.85,6,20,0); 
filtered_stimTrial_paddings1=filtfilthd(Hd,stimTrial_paddings);
close;figure;
for i=1:16
    plot(filtered_stimTrial_paddings1(:,:,i));
    pause(0.5);
end


% low pass 
% (this will create short edge artifact )
Hd = makefilter(sr,8,8.5,6,20,0);  % Hd = makefilter(sr,50,51,6,20,1); 
filtered_stimTrial_paddings2=filtfilthd(Hd,filtered_stimTrial_paddings1);% after calling EEGLAB filtfilt reset to matlab default
% open /home/zhibin/Documents/GitHub/matlab-archive/hnlcode/external/fieldtrip-20141007/external/signal/filtfilt.m % this not work
% /usr/local/MATLAB/R2019a/toolbox/signal/signal/filtfilt.m % this works
close;figure;
for i=1:16
    plot(filtered_stimTrial_paddings2(:,:,i));
    ylim([-50 50]);
    pause(0.5);
end

% remove paddings
filtered_stimTrial_data=filtered_stimTrial_paddings2((size(padding,1)+1):(size(padding,1)+size(stimTrial_data,1)),:,:);
clear filtered_stimTrial_paddings1 stimTrial_paddings 
close;figure
for i=1:16
    plot(filtered_stimTrial_data(:,:,i));
    ylim([-50 50]);
    pause(0.5);
end


%%
% remove 3 prefrontal channels
tapTrial_erp=mean(tapTrial_data(:,4:32,:),3);
% tapTrial_erp=mean(tapTrial_data(:,1:32,:),3);

[U,S,V]=svd(tapTrial_erp,"econ");
figure;image(S);
% image(U);image(V);
diagS=diag(S) % the eigen values (a set of weights) for the 32 channels 

% See if U(:,1) looks like the readiness potential?
% plot(-1*U(:,1));xlabel('time');title('U(:,1)');xline(2000,'m',{'reponse'});
for i=1:4
    figure;
    plot(-1*U(:,i));xlabel('time');
    title(['-1*U(:,' num2str(i) '); eigen value ' num2str(i) ' = ' num2str(diagS(i))]);xline(2000,'m',{'reponse'});
end
% figure;plot(-1*U(:,2));xlabel('time');title('-1*U(:,2)');xline(2000,'m',{'reponse'});
% figure;plot(-1*U(:,3));xlabel('time');title('-1*U(:,3)');xline(2000,'m',{'reponse'});
% figure;plot(-1*U(:,4));xlabel('time');title('-1*U(:,4)');xline(2000,'m',{'reponse'});

% topoplot to examine them
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/channels_info
load('chaninfo.mat')
% chaninfo29=chaninfo(4:32);
figure;
% topoplot(V(:,1),chaninfo,'nosedir','+X');
% topoplot(V(:,1),chaninfo29,'nosedir','+X');title('+1*V(:,1)');
for i=1:2
    figure;
    topoplot(-1*V(:,i),chaninfo29,'nosedir','+X');title(['-1*V(:,' num2str(i) ')']);
end

% stimTrials_synchronization_erp=mean(stimTrials_synchronization_r,3);
stimTrial_erp_estimate_trial=zeros(size(filtered_stimTrial_data,1),1);
for i=1:240
    stimTrial_erp_estimate_trial(:,i)=squeeze(filtered_stimTrial_data(:,4:32,i))*(-1*V(:,1));
end

% % Examine the single trial estimate (take time)
% close; figure('units','normalized','outerposition',[0 0 1 1]);
% for i=1:240
%     subplot(15,16,i)
%     plot(stimTrial_erp_estimate_trial(:,i));
%     title(['trial' num2str(i)]); xline(2000,'m',{'task'});
% end
% suptitle('single trial estimate of ERP (multiply with -1*V(:,1), look for negative peaks)');

% Examine 16 of them
close; figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:16
    subplot(4,4,i)
    plot(stimTrial_erp_estimate_trial(:,i));
    title(['trial' num2str(i)]); xline(2000,'m',{'task'});
end
suptitle('single trial estimate of ERP (multiply with -1*V(:,1), look for negative peaks)');

% % Plot together (hard to see)
% close; figure('units','normalized','outerposition',[0 0 1 1]);
% plotx(stimTrial_erp_estimate_trial(:,1:16));
% xline(2000,'m',{'task'});
% title('single trial estimate of ERP (multiply with -1*V(:,1), look for negative peaks)');
% ylim([-20 20])
%% Find the negative minima and save the EEG errors
minimal_inds=zeros(1,240);
error=zeros(1,240);
for i=1:240
    minimal_value=min(stimTrial_erp_estimate_trial(:,i));
    minimal_ind=find(stimTrial_erp_estimate_trial(:,i)==minimal_value);
    minimal_inds(i)=minimal_ind;
    error(i)=minimal_ind-2000;
end
plot(minimal_inds,'ro');
plot(error,'ro');

%% mark and plot the errors
% Examine 16 of them
close; figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:16
    subplot(4,4,i)
    plot(stimTrial_erp_estimate_trial(:,i));
    title(['trial' num2str(i)]); xline(2000,'m',{'task'});
    hold on; xline(minimal_inds(i),'k',{'minima'});
end
suptitle('single trial estimate of ERP < 8 Hz(x -1*V(:,1), look for negative peaks)');

% % Examine the single trial estimate (take time)
% close; figure('units','normalized','outerposition',[0 0 1 1]);
% for i=1:240
%     subplot(15,16,i)
%     plot(stimTrial_erp_estimate_trial(:,i));
%     title(['trial' num2str(i)]); 
%     xline(2000,'m');
%     xline(minimal_inds(i),'k');
% end
% suptitle('single trial estimate of ERP (multiply with -1*V(:,1), look for negative peaks)');

%% 3 plots
y=error/2000;
% % remove error > 1 s
% indices=find(y>1);
% y(indices)=NaN;
% y(indices)=[];
% % remove error < -1 s
% indices=find(y<-1);
% y(indices)=NaN;
% y(indices)=[];
% % Error1(indices)=2;

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
    
figure('units','normalized','outerposition',[0 0 1 0.3]);

subplot(1,3,1);
plot(y,'k.');
xlabel('Trial');ylabel('Error(sec.)');title('A. Time Series');

subplot(1,3,2);
autocorr(y,'NumLags',50,'NumSTD',2);
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
    
