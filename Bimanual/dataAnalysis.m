% the percentage of time in the square 
% PercentIn=sum(NumInside)/length(Pos); 

% plot the mouse trace
for t=1:5
    xLyLt=dataL(t).xLyL;
    xRyRt=dataR(t).xRyR;
    
    figure;
    
    plot(xLyLt(:,1),xLyLt(:,2),'-b');
    set(gca, 'YDir', 'reverse');
    hold on;
    plot(PosL(1,:),PosL(2,:),'-r'); 
    
    xlim([0 screenXpixels]);ylim([0 screenYpixels]);
    title(['trial' num2str(t)]);
    
    plot(xRyRt(:,1),xRyRt(:,2),'-b');
    set(gca, 'YDir', 'reverse');
    hold on;
    plot(PosR(1,:),PosR(2,:),'-r'); 
    
end