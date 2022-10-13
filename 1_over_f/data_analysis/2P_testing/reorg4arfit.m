function [Int12LR] = reorg4arfit(test_data)
% Reorganized the data for arfit function input in a Int12LR matrix of
% samples x 2LR x ntr

% Method 1: truncate the time series to 100 samples
% Int12LR=nan(100,2,size(test_data,2)*size(test_data,3));
% nL=1;nR=1;
% for p=1:size(test_data,1)
%     for b=1:size(test_data,2)
%         for s=1:size(test_data,3)
%             if p==1;
%                 Int12LR(1:100,1,nL)=test_data{1,b,s}(1:100,:);
%                 nL=nL+1;
%             end
%             if p==2
%                 Int12LR(1:100,2,nR)=test_data{2,b,s}(1:100,:);
%                 nR=nR+1;
%             end
%         end
%     end
% end


% resample to 200 samples
Maxlength=200;
Int12LR=nan(Maxlength,2,size(test_data,2)*size(test_data,3));
nL=1;nR=1;
for p=1:size(test_data,1)
    for b=1:size(test_data,2)
        for s=1:size(test_data,3)
            if p==1;
                Int12LR(1:Maxlength,1,nL)=resample(test_data{1,b,s},Maxlength,length(test_data{1,b,s}));
                nL=nL+1;
            end
            if p==2
                Int12LR(1:Maxlength,2,nR)=resample(test_data{2,b,s},Maxlength,length(test_data{2,b,s}));
                nR=nR+1;
            end
        end
    end
end


end