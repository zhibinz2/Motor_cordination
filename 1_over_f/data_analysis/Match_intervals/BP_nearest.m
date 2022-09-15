function [samples, intervals, BPL, BPR, intL, intR] = BP_nearest(BPL,BPR, tooshort)


sampL = find(BPL==1);
sampR = find(BPR==1);
intL=diff(sampL);
intR=diff(sampR);
nL = length(sampL);
nR = length(sampR);
doubleL = find(intL < tooshort);
doubleR = find(intR < tooshort);
if ~isempty(doubleL) 
	BPL(sampL(doubleL+1)) = -1;
end
if ~isempty(doubleR)
	BPR(sampL(doubleR+1)) = -1;
end
sampL = find(BPL>0);
sampR = find(BPR>0);
intL=diff(sampL);
intR=diff(sampR);
nL = length(sampL);
nR = length(sampR);
% nR is shorter this time
dist = [];
for j = 1:length(sampL) 
	dist(j,:) = sampL(j) - sampR;
end;

[mindistL,mindistindexL] = min(abs(dist),[],2);% corresponding to each tap from L
% mindistindexL is the index of the tap with minimal distance in samp R (same length as L)
[mindistR,mindistindexR] = min(abs(dist),[],1);% corresponding to each tap from R
% mindistindexR is the index of the tap with minimal distance in samp L (same length as R)

%check for redundancy 
ntapL = length(mindistL);
ntapR = length(mindistR);

if ntapL <= ntapR 
	badlist = zeros(1,ntapL);
	[tapuniqueL,ia,ic] = (unique(mindistindexL));
	if length(tapuniqueL) ~= ntapL 
		repeats = find(diff(ic) == 0);
		diffrep = mindistL(repeats)-mindistL(repeats+1);
		repeats(diffrep<0) = repeats(diffrep<0)+1;		
		badlist(repeats) = 1; 
	end;
	goodlist = setdiff([1:ntapL],find(badlist)); % for L % or ~badlist 
	goodlistR = mindistindexL(goodlist);
	if goodlist(1) == 1 | goodlistR(1) == 1
		goodlistR = goodlistR(2:end);
		goodlist = goodlist(2:end);
	end
	samples(:,1) = sampL(goodlist);
	samples(:,2) = sampR(goodlistR);
	intervals(:,1) = intL(goodlist-1);
	intervals(:,2) = intR(goodlistR-1);
elseif ntapR < ntapL % will remove bad tap in R
	badlist = zeros(1,ntapR); % will mark tap to be removed as 1 in R
	[tapuniqueR,ia,ic] = (unique(mindistindexR));
    % ic has the same length as mindistindexR. (same length as R)
    % minidistindexR=tapuniqueR(ic).
	if length(tapuniqueR) ~= ntapR % meaning at least one same tap in L has minimal distance with multiple taps from R
		repeats = find(diff(ic) == 0); % index of the 1st repeat in R with minial dist to the same tap in L
		diffrep = mindistR(repeats)-mindistR(repeats+1); % this is between 2 local repeats everyhwere
		repeats(diffrep<0) = repeats(diffrep<0)+1; % diffrep<0 means the minimal distance increased in the next tap of R	 
		badlist(repeats) = 1; % if minimal distance increased in the next tap of R, mark this next tap as bad
	end;
	% goodlist = setdiff([1:ntapR],badlist);% Use this to extract interval as intR(goodlist)
	goodlist = setdiff([1:ntapR],find(badlist));% The above line is wrong? %or  ~badlist
    goodlistL = mindistindexR(goodlist); 
	if goodlist(1) == 1 | goodlistL(1) == 1
		goodlistL = goodlistL(2:end);
		goodlist = goodlist(2:end);
	end;
	samples(:,1) = sampL(goodlistL);% Use this for EEG extraction
	samples(:,2) = sampR(goodlist);
	intervals(:,1) = intL(goodlistL-1); % Use this 
	intervals(:,2) = intR(goodlist-1); 
end;

