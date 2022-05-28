%% resample key presses
BPCondi1L;BPCondi1R;
figure;plot(Calinterval(BPCondi1L'),'b.');hold on;
Addshowframes=Calinterval(BPCondi1L');
Addshowframes=cumsum(Addshowframes);
plot(Addshowframes,'mo');hold on;
Addshowframes=round(Addshowframes/(2000/100)/2);
plot(Addshowframes,'m+');hold on;
Addshowframes(1) 
NumFramesInterval2Hz*29+1
shift=Addshowframes(1)-(NumFramesInterval2Hz*29+1)
Addshowframes=Addshowframes-shift;
Addshowframes(1)
NumFramesInterval2Hz*29+1
Showframes1=[];
Showframes1=[1:NumFramesInterval2Hz:NumFramesInterval2Hz*30 Addshowframes(2:end)]; 
plot(Showframes1,'go');hold on;
Showframes1=Showframes1*1.5;
plot(Showframes1,'k*');
%% Randomization of the pace stimuli in each condition ****************************************
    % Mean stimulus interval for 2Hz pacing
    MeanTapInterval2Hz=0.5; % second
    NumFramesInterval2Hz=round(MeanTapInterval2Hz/(ifi*waitframes));  
    
    % condition 1-4 (paced the frist 30 taps, set the rest of the frames with value zeros)
    Showframes1=[1:NumFramesInterval2Hz:NumFramesInterval2Hz*30 zeros(1,570)]; 
    Showframes2=Showframes1;
    Showframes3=Showframes1;
    Showframes4=Showframes1;
    
    % Mean stimulus interval for 3Hz pacing
    MeanTapInterval3Hz=1/3; % second
    NumFramesInterval3Hz=round(MeanTapInterval3Hz/(ifi*waitframes));  % on average 72 frames per stimulus 
    % condition 5
    Showframes5=[1:NumFramesInterval3Hz:NumFramesInterval3Hz*30 zeros(1,570)]; 
    Showframes6=Showframes5;
    
    % combine all 6 conditions
    Showframes=[Showframes1;Showframes2;Showframes3;Showframes4;Showframes5;Showframes6];