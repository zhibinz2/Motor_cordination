% This script only works on ALIENWARE PC at home

cd D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\EEGanalysis\20211102

detrend_data=ndetrend(transpose_data,1); 

sgolayfilt(detrend_data,1,50);