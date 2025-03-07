function [samples, intervals, BPL, BPR, intL, intR] = BP_nearest(BPL,BPR, tooshort)
sampL = find(BPL==1); % sample index of the tap
sampR = find(BPR==1);
intL=diff(sampL); % intervals
intR=diff(sampR);
nL = length(sampL); % number of taps
nR = length(sampR);
doubleL = find(intL < tooshort); % index of bad tap (double tapping)
doubleR = find(intR < tooshort);
if ~isempty(doubleL) 
	BPL(sampL(doubleL+1)) = -1; % mark double tapping as -1 in the sample index
end
if ~isempty(doubleR)
	BPR(sampR(doubleR+1)) = -1; % change sampL to sampR on 10/04/2022 AM possibly need to redo all previous data?
end
sampL = find(BPL>0); % keep all good taps so far
sampR = find(BPR>0);
intL=diff(sampL); % update the intervals
intR=diff(sampR);
nL = length(sampL); % update the number of taps
nR = length(sampR);
% nR is shorter this time
dist = []; % distance between each tap from L and R in a matrix of nL X nR
for j = 1:length(sampL) 
	dist(j,:) = sampL(j) - sampR;
end

% find the nearest tap from the other 
[mindistL,mindistindexL] = min(abs(dist),[],2);% corresponding to each tap from L
% mindistL is nearest distance list to each tap in L
% mindistindexL is the index of the tap with minimal distance to samp R (same length as R)
[mindistR,mindistindexR] = min(abs(dist),[],1);% corresponding to each tap from R
% mindistR is nearest distance list to each tap in R
% mindistindexR is the index of the tap with minimal distance to samp L (same length as L)


% check for redundancy ###################################################
ntapL = length(mindistL); % the number of tap distances (same length as R)
ntapR = length(mindistR);

if ntapL <= ntapR % This means nL>nR. Will remove bad tap in L
	badlist = zeros(1,ntapL); % initialize badlist (will mark tap to be removed in L as 1)
	[tapuniqueL,ia,ic] = (unique(mindistindexL)); % find unique values (remove redundancy)
    % ic has the same length as mindistindexL. (same length as R)
    % minidistindexL=tapuniqueL(ic).
    % mindistindexL(ia)=tapuniqueL. ia is not used
	if length(tapuniqueL) ~= ntapL % then there is redundendy in mindistindexL, meaning at least one same tap in R has minimal distance with multiple taps from L
		repeats = find(diff(ic) == 0);  % index of the 1st repeat in L with minial dist to the same tap in R
		diffrep = mindistL(repeats)-mindistL(repeats+1); % this is between the 2 local repeats everyhwere
		repeats(diffrep<0) = repeats(diffrep<0)+1; % diffrep<0 means the minimal distance increased in the next tap of L
		badlist(repeats) = 1; % if minimal distance increased in the next tap of L, mark this next tap as bad
    end
	goodlist = setdiff([1:ntapL],find(badlist)); % for L % or ~badlist 
	goodlistR = mindistindexL(goodlist); % get the correponding index in R
	if goodlist(1) == 1 | goodlistR(1) == 1 % if the first tap is good
		goodlistR = goodlistR(2:end); % Then start list from the second tap, because index of the intervals are counted before each tap from the second tap.
		goodlist = goodlist(2:end);
	end
	samples(:,1) = sampL(goodlist);
	samples(:,2) = sampR(goodlistR);
	intervals(:,1) = intL(goodlist-1);
	intervals(:,2) = intR(goodlistR-1);
elseif ntapR < ntapL % will remove bad tap in R
	badlist = zeros(1,ntapR); % will mark tap to be removed as 1 in R
	[tapuniqueR,ia,ic] = (unique(mindistindexR)); % find unique values (remove redundancy)
    % ic has the same length as mindistindexR. (same length as L)
    % minidistindexR=tapuniqueR(ic).
    % mindistindexR(ia)=tapuniqueR. ia is not used
	if length(tapuniqueR) ~= ntapR % meaning at least one same tap in L has minimal distance with multiple taps from R
		repeats = find(diff(ic) == 0); % index of the 1st repeat in R with minial dist to the same tap in L
		diffrep = mindistR(repeats)-mindistR(repeats+1); % this is between 2 local repeats everyhwere
		repeats(diffrep<0) = repeats(diffrep<0)+1; % diffrep<0 means the minimal distance increased in the next tap of R	 
		badlist(repeats) = 1; % if minimal distance increased in the next tap of R, mark this next tap as bad
    end
	% goodlist = setdiff([1:ntapR],badlist);% Use this to extract interval as intR(goodlist)
	goodlist = setdiff([1:ntapR],find(badlist));% The above line is wrong? %or  ~badlist
    goodlistL = mindistindexR(goodlist); 
	if goodlist(1) == 1 | goodlistL(1) == 1
		goodlistL = goodlistL(2:end);
		goodlist = goodlist(2:end);
    end
	samples(:,1) = sampL(goodlistL);% Use this for EEG extraction
	samples(:,2) = sampR(goodlist);
	intervals(:,1) = intL(goodlistL-1); % Use this 
	intervals(:,2) = intR(goodlist-1); 
end

