%% Plot an illustration of the photocells

%     event types   
%     16384 (upper left)    
%     32768 (upper right)     
%     1200001 (upper right)    
%     1200002 (upper left)

%% Photocells illustration during Resting state EEG acquisition
% Upper left and bottom right (4 times)

close;figure;
subplot(2,2,1);
plot([0 1.5 2 3.5],[1 1 1 1],'ro');xlim([-0.5 4]);ylim([0 2]);
xlabel('time(min)');ylabel('Upper Left Photocell');
xline(0,'g-','open eye start');xline(1.5,'g-','open eye start');
xline(2,'k-','closed eye start');xline(3.5,'k-','closed eye start');
yline(1,'m--','time line');
title('(event type 1200002 or 16384)');
subplot(2,2,2);
plot([1000 1500 2500],[NaN NaN NaN],'ro');xlim([-0.5 4]);ylim([0 2]);
xlabel('time(min)');ylabel('Upper Right Photocell');
xline(0,'g-','open eye start');xline(1.5,'g-','open eye start');
xline(2,'k-','closed eye start');xline(3.5,'k-','closed eye start');
yline(1,'m--','time line');
title('(event type 1200001 or 32768)');
subplot(2,2,3);
plot(linspace(1000,2500,45),NaN(1,45),'ro');xlim([-0.5 4]);ylim([0 2]);
xline(0,'g-','open eye start');xline(1.5,'g-','open eye start');
xline(2,'k-','closed eye start');xline(3.5,'k-','closed eye start');
yline(1,'m--','time line');
xlabel('time(min)');ylabel('Bottom Left Photocell');
title('peaks in analog1=133');
subplot(2,2,4);
plot([0 1.5 2 3.5],[1 1 1 1],'ro');xlim([-0.5 4]);ylim([0 2]);
xline(0,'g-','open eye start');xline(1.5,'g-','open eye start');
xline(2,'k-','closed eye start');xline(3.5,'k-','closed eye start');
yline(1,'m--','time line');
xlabel('time(min)');ylabel('Bottom Right Photocell');
title('peaks in analog2=134');


%% photocells illustration during trials
close;figure;
subplot(2,2,1);
plot([1000 1500 2500],[1 NaN 1],'ro');xlim([0 3000]);ylim([0 2]);
xlabel('time(ms)');ylabel('Upper Left Photocell');
xline(1000,'g-','trial start');xline(1500,'g-','movement start');xline(2500,'k-','trial end');
yline(1,'m--','time line');
title('(event type 1200002 or 16384)');
subplot(2,2,2);
plot([1000 1500 2500],[NaN 1 NaN],'ro');xlim([0 3000]);ylim([0 2]);
xlabel('time(ms)');ylabel('Upper Right Photocell');
xline(1000,'g-','trial start');xline(1500,'g-','movement start');xline(2500,'k-','trial end');
yline(1,'m--','time line');
title('(event type 1200001 or 32768)');
subplot(2,2,3);
plot(linspace(1000,2500,45),ones(1,45),'ro');xlim([0 3000]);ylim([0 2]);
xline(1000,'g-','trial start');xline(1500,'g-','movement start');xline(2500,'k-','trial end');
yline(1,'m--','time line');
xlabel('time(ms)');ylabel('Bottom Left Photocell');
title('peaks in analog1=133');
subplot(2,2,4);
plot([1000 1500 2500],[1 1 1],'ro');xlim([0 3000]);ylim([0 2]);
xline(1000,'g-','trial start');xline(1500,'g-','movement start');xline(2500,'k-','trial end');
yline(1,'m--','time line');
xlabel('time(ms)');ylabel('Bottom Right Photocell');
title('peaks in analog2=134');