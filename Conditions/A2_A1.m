function [xL,yL,CentroidL,PosL,PosL0,PosL1,PosL2,PosL3,ThicknessL,ColorL,xR,yR,CentroidR,PosR,PosR0,PosR1,PosR2,PosR3,ThicknessR,ColorR]=A2_A1(steplength,yCenter,screenXpixels)

%UNTITLED8 Summary of this function goes here
[xL,yL]=drawA2(steplength); 
[xR,yR]=drawA1(steplength);
[CentroidL,PosL,PosL0,PosL1,PosL2,PosL3,ThicknessL,ColorL] = DrawPosL(xL,yL,steplength,yCenter);
[CentroidR,PosR,PosR0,PosR1,PosR2,PosR3,ThicknessR,ColorR] = DrawPosR(xR,yR,steplength,yCenter,screenXpixels);

end

