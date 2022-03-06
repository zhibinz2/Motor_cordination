

%% https://www.fieldtriptoolbox.org/faq/should_i_add_fieldtrip_with_all_subdirectories_to_my_matlab_path/
addpath C:\Users\zhibi\Downloads\fieldtrip-lite-20220304\fieldtrip-20220304
ft_defaults
fieldtripdefs.m


%%  https://www.fieldtriptoolbox.org/tutorial/nirs_singlechannel/
help ft_preprocessing
edit ft_preprocessing

cfg = [];
cfg.dataset = 'C:\Users\zhibi\Downloads\nirs_singlechannel\motor_cortex.oxy3';

[data] = ft_preprocessing(cfg);


%% https://www.fieldtriptoolbox.org/tutorial/nirs_multichannel/