function [BP_L] = cutofflast(BP_L)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
ind=find(BP_L==1);
BP_L=BP_L(1:ind(end));
end

