function inside = IsInDot(x,y,Position,ConnectDotSize)
% 
% Is location x,y inside the circle/dot?
% Position is the center of the circle
% ConnectDotSize is the radius of the circle

% 20210726  zhibin wrote it

if sqrt((x-Position(1)).^2+(y-Position(2)).^2)<ConnectDotSize
	inside = 1;
else
	inside = 0;
end
