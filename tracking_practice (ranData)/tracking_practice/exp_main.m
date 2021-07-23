function exp_main

close all;
clc;
clear;

addpath data
addpath Track_core
addpath GUI

rng('shuffle');


scaleval=.1;
 
% gui stuff goes here
h = gui_init();
subinfo.subcode = get(h.hEdt1,'String');
subinfo.subage = str2double(get(h.hEdt2,'String'));
subinfo.subear = get(h.hEdt3,'Value');
subcode = subinfo.subcode;
[protocol, Subdata] = load_subject_info(subcode);
if protocol ==9999
    disp('Data collection is done for this subject.');
    return;
end
% chn = subinfo.subear;   %this info7rmation is ignored

set(h.hLbl4,'String',{['Condition ' num2str(protocol(1))]});

    
stop_flag = 0;
while stop_flag ==0

dur=0.6;  fs=44100; %duration in seconds; 300 ms is stimulus
    ntrials = 50;
    track = [];
    [x, tracking_par] = Track_init();    %x is deltaf/f percent
    
    for itrial = 1:ntrials
        [answ,resp,cor,cuefund,tonefund] = run_trial(x,h,protocol,scaleval,itrial);
        fdbk(cor, h.hText);
        tracking_par.cor = cor;
        [dx, rev_flag, tracking_par] = Track_update(tracking_par);
               
        track(itrial,:) = [x, answ, resp, cor, rev_flag,cuefund,tonefund];
        x = max(min(x + dx,max(tracking_par.tracklimit)),min(tracking_par.tracklimit));
        %             x = x + dx;
        pc = round(sum(track(1:itrial,4))/itrial*100);
        set(h.hLbl4,'String',{['Condition ' num2str(protocol(1)) ...
            ', Trial: ' num2str(itrial) '/' num2str(ntrials)] ;...
            ['PC=' num2str(pc) ' percent']});
    end
     Subdata.data{end+1}.track = track;
     Subdata.data{end}.param = tracking_par;
     Subdata.data{end}.subinfo = subinfo;
     Subdata.data{end}.cond = protocol(1);
     Subdata.data{end}.timestamp = now;
    

    save([cd '\data\', subcode, 'data.mat'], '-struct', 'Subdata');


    set(h.hText,'String','Press "q" to quit, press any other key to continue.');
    while waitforbuttonpress==0
    end
    set(h.hText,'String','Please wait');
    drawnow;
    charval = get(h.hFig,'CurrentCharacter');
    if strcmp(charval,'q')
        stop_flag =1;
        close(h.hFig);
    else
        pause(3);
    end
    protocol = load_subject_info(subcode);  % load conditions for the next track
    if protocol ==9999
        stop_flag = 1;
    end
    
end

set(h.hText,'String','Finished!');

end
