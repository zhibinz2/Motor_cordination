%loop_thresh.m   
% quick program to estiamte thresholds for differnet conditions
clear all 
close all

addpath ../data

clc

disp(' ');
disp(' to estimate thresholds ');
obs=input('Enter observer code ','s');
strval=[obs 'data'];
load(strval);
clear strval


for iblk=1:length(data)   
    data_tmp = data{iblk};
    cond(iblk)=data_tmp.cond;
    track = data_tmp.track;
    rev_idx = find(track(:,5)==1);
    if length(rev_idx) < 6
        thre_tmp(iblk) = 999;
    else
        thre_tmp(iblk) = mean(track(rev_idx(end-5:end),1));
    end
    figure(1)
    plot(track(:,1))
    xlabel(' trial number '); ylevel=('delta f0/f0');
    title([obs]);
    pause(1);
end


disp('   '); disp('   '); 
disp([' observer ' obs])
disp('thresh %delfo/fo');
for iblk=1:length(data)
    disp(['      ' num2str((thre_tmp(iblk))*100) '%'])
end


