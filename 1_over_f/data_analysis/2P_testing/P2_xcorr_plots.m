%% P2 behaviral (short-range / local statistics; weak anticipation)
% [ autocorr(BP_L) + xcorr(BP_L vs FB_L) + xcorr(BP_R vs FB_R) + autocorr(BP_R) ]
% 6 conditions -subplots(4,6,i)


%% P2 behaviral (long-range / global statistics; strong anticipation)
% [ Pspectra(BP_L) + Pspectra(BP_R) ] + [DFA(BP_L) + DFA(BP_R)] 
% 6 conditions -subplots(2,6,i)


%% P2 Error (short & long -range statistics)
% autocorr (Error) + Pspectra(Error) + DFA (Error)
% 6 conditions -subplots(3,6,i)


% use a sliding window and do it again
% 6 conditions -subplots(3,6,i)


%% P2 EMG (short-range / local statistics; weak anticipation)
% [ autocorr(EMG_L) + xcorr (EMG_L vs EMG_R) + autocorr(EMG_R) ]
% 6 conditions - subplots(3,6,i)


%% P2 EMG (long-range / global statistics; strong anticipation)
% [ Pspectra(EMG_L) + Xspectra (EMG_L vs EMG_R) + Pspectra(EMG_R) ]
% 6 conditions - subplots(3,6,i)


%% EEG auto xcorr (weak anticipation, short-range / local statistics)
% [ autocorr(EEG_L) + xcorr(EEG_L vs EEG_R) + autocorr(EEG_R) ]
% each of the 6 conditions (subplots(3,6,i)) for each of the 32 channels;
% auto save the 32 figures in a folder named -


%% EEG Power spectra and Xspectra (strong anticipation, long-range / global statistics)
% [ Pspectra(EEG_L) + Xspectra(EEG_L vs EEG_R) + Pspectra(EEG_R) ]
% each of the 6 conditions (subplots(3,6,i)) for each of the 32 channels;
% auto save the 32 figures in a folder named - 





