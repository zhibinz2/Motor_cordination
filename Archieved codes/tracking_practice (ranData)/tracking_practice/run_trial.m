function [answ,resp,cor,cuefund,tonefund] = run_trial(x,h,protocol,scaleval,itrial)

% runs tracking for cue-tone

fs = 44100;
nAFC = 2;
answ = randi(nAFC);   % 1 (down) or 2(up)
dur=0.6;   %duration in seconds; 250 ms is stimulus
isi=0.2;  % will pad the stimulus 
gatedur=0.01;  %gate in seconds

t=0:1/fs:(dur-(1/fs));
pauser=zeros(1,6000);

% gate
whichTime = t <= gatedur;
gateTimes = t(whichTime);
rise = (1-cos(pi*gateTimes/gatedur))/2;
fall = flip(rise);
tails =  length(rise);
middle = ones(1,(length(t) - 2*tails));
gate = [rise middle fall];

% choosing a cue/tone fundfreq between 600-1000 Hz
while 1
    a = rand(1,1)*(log10(1000)-log10(600))+log10(600);
    cuefund = 10^a;
    
    % tone fundfreq
    if answ == 1 % answer is "low"
        tonefund = (1-x)*cuefund;
    elseif answ == 2 % answer is "high"
        tonefund = (1+x)*cuefund;
    end
    
    % tonefund should be between 600 and 1000
    if tonefund > 600 && tonefund < 1000
        break
    end
end

% missing fundamental harmonics
lowHarmo = randi([2 4],1); % lowest harmonic number in the cue, 2, 3 or 4
harmonics = zeros(1,fs*dur);
harmophase = 2*pi*rand(1,3);
k = 1;
for i = lowHarmo:lowHarmo+2
    harmonics = harmonics + scaleval*cos(2*pi*(cuefund*i)*t+harmophase(k)); % cue
    k = k+1;
end
harmonics = [pauser gate.*harmonics];

% tone
tone = zeros(size(t));
tonephase = 2*pi*rand(1,5);
if protocol(1) == 1
    tone = scaleval*cos(2*pi*tonefund*t+tonephase(3));
end

tone = [pauser gate.*tone];

player = audioplayer(harmonics, fs,24);
pause(0.1)
set(h.hText,'string','Cue');
playblocking(player);
pause(0.01)    
set(h.hText,'string',[]);
pause(isi);   %most of ISI

player = audioplayer(tone, fs,24);
pause(0.1)
set(h.hText,'string','Tone');
playblocking(player);
pause(0.1);  
set(h.hText,'string',[]);

pause(0.01);
set(h.hText,'string',...
    {'Please select the Hi-Lo (1) or Lo-Hi (2) ',...
    'Press "1" or "2".'});

resp = 0;
while all(resp~=[1 2])
    while waitforbuttonpress==0
    end
    set(h.hText,'string',' ');
    charval = get(h.hFig,'CurrentCharacter');
    resp = str2double(charval);
end

if resp == answ
    cor = 1;
else
    cor = 0;
end