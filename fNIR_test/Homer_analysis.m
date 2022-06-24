%% Homer3
cd /home/zhibin/Documents/GitHub/Motor_cordination/fNIR_test

cd /home/zhibin/Downloads/homer3_install_v1_33_0_src/Homer3/
run setpaths.m
cd /home/zhibin/Documents/Artinis_NIRS/shane/20220623
addpath(genpath('/home/zhibin/Documents/GitHub/oxysoft2matlab_v1.81'));
oxysoft2matlab % select both 20220507.oxy3 and 20220507.oxyproj files
% convert oxy5 to snirf
Homer3


%% Homer2
cd /home/zhibin/Downloads/homer2_src_v2_8_11022018/homer2/
run setpaths.m