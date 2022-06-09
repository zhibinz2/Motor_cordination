ifi
1/ifi
mean(fliptimes)
max(fliptimes)
min(fliptimes)
%%
close all;
figure('units','normalized','outerposition',[0 0 0.4 0.5]);
plot(fliptimes*1000,'k*'); 
title(['Flip once every ' num2str(waitframes) ' video frame (screen refresh rate: ' num2str(round(1/ifi)) ' Hz)'])
yline(waitframes*ifi*1000,'m',{['Time to flip every ' num2str(waitframes) ' video frames']})
ylabel(['Actual flip interval: ' num2str(floor((min(fliptimes))*1000)) ' - ' num2str(ceil((max(fliptimes))*1000)) ' ms'])
xlabel('Flips')
ylim([0 max(fliptimes)*1000+5])

%% Left Player
% Left player  (read button press)
% initialize the time
pressedL1=[];RBkeyL1=[];
Start_times=now;
while pressedL1 ~= 1
        [pressedL1, ~]=readCedrusRB(deviceL, keymapL) % extract first key press
end
pressedL1
End_time=now;
Get_keypress_time=(End_time-Start_times)* 10 ^ 5
num2str(Get_keypress_time)
