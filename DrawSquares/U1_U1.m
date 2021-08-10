function [xL,yL,PosL,PosL0,PosL1,PosL2,PosL3,ThicknessL,ColorL,xR,yR,PosR,PosR0,PosR1,PosR2,PosR3,ThicknessR,ColorR]=U1_U1(steplength,yCenter,screenXpixels)

%UNTITLED8 Summary of this function goes here
[xL,yL]=drawU1(steplength); 
[xR,yR]=drawU1(steplength);
[PosL,PosL0,PosL1,PosL2,PosL3,ThicknessL,ColorL] = DrawPosL(xL,yL,steplength,yCenter);
[PosR,PosR0,PosR1,PosR2,PosR3,ThicknessR,ColorR] = DrawPosR(xR,yR,steplength,yCenter,screenXpixels);

end

