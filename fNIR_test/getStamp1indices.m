function [stamp1index_Baseline_Start,stamp1index_Baseline_End,stamp1index_Trial_Start,stamp1index_Trial_End] = getStamp1indices(time_stamps1,time_stamps2,time_series2);
% indexing the Markers in timestamp1 base on markers in timestamp2
%   Detailed explanation goes here
index_Baseline_Start = find(strcmp(horzcat({time_series2{:}}),'Baseline-Start')==1);
index_Baseline_End = find(strcmp(horzcat({time_series2{:}}),'Baseline-End')==1);
index_Trial_Start = find(strcmp(horzcat({time_series2{:}}),'Trial-Start')==1);
index_Trial_End = find(strcmp(horzcat({time_series2{:}}),'Trial-End')==1);
index_Screen_Flip = find(strcmp(horzcat({time_series2{:}}),'Screen_Flip')==1);

repmat_stamp1_Baseline_Start=repmat(time_stamps1',[1,length(time_stamps2(index_Baseline_Start))]);
[minValue_Baseline_Start,stamp1index_Baseline_Start] = min(abs(repmat_stamp1_Baseline_Start-time_stamps2(index_Baseline_Start)));

repmat_stamp1_Baseline_End=repmat(time_stamps1',[1,length(time_stamps2(index_Baseline_End))]);
[minValue_Baseline_End,stamp1index_Baseline_End] = min(abs(repmat_stamp1_Baseline_End-time_stamps2(index_Baseline_End)));

repmat_stamp1_Trial_Start=repmat(time_stamps1',[1,length(time_stamps2(index_Trial_Start))]);
[minValue_Trial_Start,stamp1index_Trial_Start] = min(abs(repmat_stamp1_Trial_Start-time_stamps2(index_Trial_Start)));

repmat_stamp1_Trial_End=repmat(time_stamps1',[1,length(time_stamps2(index_Trial_End))]);
[minValue_Trial_End,stamp1index_Trial_End] = min(abs(repmat_stamp1_Trial_End-time_stamps2(index_Trial_End)));
end

