% Calibration.m
addpath analysis

load('calibration.mat');
obj = audioplayer(y,fs);
play(obj)
