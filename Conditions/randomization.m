A1_A2
A1_A1
A1_U2
A1_U1

A2_A1
A2_A2
A2_U1
A2_U2

U1_U2
U1_U1
U1_A2
U1_A1

U2_U1
U2_U2
U2_A1
U2_A2

%%
'A1_A2'
'A1_A1'
'A1_U2'
'A1_U1'

'A2_A1'
'A2_A2'
'A2_U1'
'A2_U2'

'U1_U2'
'U1_U1'
'U1_A2'
'U1_A1'

'U2_U1'
'U2_U2'
'U2_A1'
'U2_A2'

%%
addpath Bimanual

conditionfuctions = {@A1_A2, @A1_A1, @A1_U2, @A1_U1, ...
    @A2_A1, @A2_A2, @A2_U1, @A2_U2, ...
    @U1_U2, @U1_U1, @U1_A2, @U1_A1, ...
    @U2_U1, @U2_U2, @U2_A1, @U2_A2 };
conditionPermutations = randperm(16); 
for cP = conditionPermutations
    [xL,yL,PosL,PosL0,PosL1,PosL2,PosL3,ThicknessL,ColorL,...
        xR,yR,PosR,PosR0,PosR1,PosR2,PosR3,ThicknessR,ColorR] ...
        = conditionfuctions{cP}(steplength,yCenter,screenXpixels);
end


%%
A1_A2 % run A1_A2.m
[xL,yL]=drawA1(steplength); 
[xR,yR]=drawA2(steplength);
[PosL,PosL0,PosL1,PosL2,PosL3,ThicknessL,ColorL] = DrawPosL(xL,yL,steplength,yCenter);
[PosR,PosR0,PosR1,PosR2,PosR3,ThicknessR,ColorR] = DrawPosR(xR,yR,steplength,yCenter,screenXpixels);

A1_A1
run 
[xL,yL]=drawA1(steplength); [xR,yR]=drawA1(steplength);
A1_U2
A1_U1

A2_A1
A2_A2
A2_U1
A2_U2

U1_U2
U1_U1
U1_A2
U1_A1

U2_U1
U2_U2
U2_A1
U2_A2