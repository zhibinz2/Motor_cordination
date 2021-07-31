% the percentage of time in the square 
% PercentIn=sum(NumInside)/length(Pos); 

% plot the mouse trace
for t=1:3
    xyt=data(t).xy;
    figure;
    plot(xyt(:,1),xyt(:,2),'-b');
    set(gca, 'YDir', 'reverse');
    hold on;
    plot(Pos(1,:),Pos(2,:),'-r'); 
    xlim([0 screenXpixels]);ylim([0 screenYpixels]);
    title(['trial' num2str(t)]);
    close all;
end