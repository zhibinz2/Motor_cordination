% psychfunc.m
clc;clear;
% Produce a psychometric function

addpath ../data;
prot_name = input('Enter observer code ','s');
ncond = input('Enter condition number ');
temp = [];

load([prot_name 'data']);
idx = length(data);

for kk = 1:idx
    if data{kk}.cond(1) == ncond
        temp = [temp;data{kk}.track(:,[1 4])]; % [x cor]
    else
        temp = temp;
    end
end

deviation = unique(temp(:,1));
pc = [];

for mm = 1:length(deviation)
    idx2 = find(temp(:,1)==deviation(mm));
    correct = sum(temp(idx2,2));
    pc(mm) = correct/length(idx2);
end

deviation = deviation';

plot(deviation,pc,'-')
title(['Psychometric Function for ',num2str(prot_name)])
xlabel('delta f0/f0')
ylabel('proportion correct')