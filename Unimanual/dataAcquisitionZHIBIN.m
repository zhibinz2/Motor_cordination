% testZhibin_Photocell_FINAL_FINAL_FINAL_matlab

% plot analog opticfiber signal
analog1=133; analog2=134;
analog1data=data(133,:);
analog2data=data(134,:);
plot(analog1data);
plot(analog2data);
close all

Halfhigh1=1/4*(max(analog1data)-min(analog1data));
Halfhigh2=1/4*(max(analog2data)-min(analog2data));
% 
ind1=find(analog1data>Halfhigh1);
plot(analog1data);
hold on;
xline(ind1(1),'--r');
xline(ind1(end),'--r');

ind1(end)-ind1(1)
ifi*numFrames*3+ifi*30*2
length(xy)*3*ifi+ifi*30*2

ifi*numFrames



% pks = findpeaks(analog1data)

pksreal=find(pks>
plot(pks,'ro');

% plot photosensor signal from stim tracker
events(:,50:100);
diff(events(1,50:100))
events(:,end-10:end)

%%
% testZhibin_Photocell_THEFINAL_matlab

plot(data(134,:))
plot(data(133,:))

events
events(:,end-10:end)
diff(events(:,end-10:end))
events(:,end-10:end)
diff(events(1,end-10:end))