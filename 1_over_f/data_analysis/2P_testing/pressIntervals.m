A=[0 0 0 0 1 1 0 0 1 1 1 1 0 0 0 0 1 1 1 1 1 0 0 0 ];
pressInd=find([0 diff(A)]==1)
pressIntervals=[pressInd(1) diff(pressInd)]