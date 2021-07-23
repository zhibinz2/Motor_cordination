function [protocol, Subdata] = load_subject_info(subcode)

if ~exist([subcode '_protocol.txt'])
    error('The protocol for this subject was not found!');
end

protocol = importdata([subcode, '_protocol.txt']);

if ~exist(['data\', subcode 'data.mat'])
    data = [];
    save(['data\', subcode 'data.mat'], 'data');
end

Subdata = load(['data\', subcode 'data.mat']);
if length(Subdata.data)<size(protocol,1)
    protocol = protocol(length(Subdata.data)+1,:);
else
    protocol = 9999;
end

end