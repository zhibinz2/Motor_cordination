function [Centroid,Pos,Pos0,Pos1,Pos2,Pos3,Thickness,Color] = DrawPosL(x,y,steplength,yCenter)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Thickness   = 15 .* ones(1,length(x)); % thickness of the route
Color  = ones(3,length(x));
% Add an offset for center the arc
x      = x + steplength/2; 
Centroidx =steplength + steplength/2; % offset circle center

% Reverse the route
% y      = abs(y-max(y)); 
% Offset the middle of the arc to be center on the middle of the screen
y      = y + steplength/2 % + yCenter - steplength/2;
Centroidy =steplength + steplength/2 % + yCenter - steplength/2; % offset circle center

% the route
Pos    = [x; y];
Centroid=[Centroidx; Centroidy];

% Postion of connecting points
Pos0   = Pos(:,1);
Pos1   = Pos(:,steplength);
Pos2   = Pos(:,steplength*2-1);
Pos3   = Pos(:,steplength*3-2);
end


%% find the center of the circle
% CenterXL = steplength + steplength/2; 
% CenterYL = steplength + (yCenter - max(y) + (max(y) - min(y))/2); 
