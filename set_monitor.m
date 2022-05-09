%% to Find Info About Your Monitors
!xrandr -q

%% to Set Your Monitor Resolution
!xrandr --output DVI-D-0 --mode 1920x1080

%% to Set Your Monitor Refresh Rate
!xrandr --output DVI-D-0 --mode 1920x1080 --rate 60.00
% for BenQ in SBSG1507 (hnlstim1)
!xrandr --output DVI-D-1 --mode 1920x1080 --rate 144.00 
!xrandr --output DP-1 --mode 1920x1080 --rate 144.00
!xrandr --output DVI-D-1 --mode 1920x1080 --rate 144.00 --output DP-1 --mode 1920x1080 --rate 144.00
% set BenQ to 60Hz  (hnlstim1)
!xrandr --output DP-1 --mode 1920x1080 --rate 60.00 --output DVI-D-1 --mode 1920x1080 --rate 60.00
% for BenQ in SBSG1503 (hnlstim2)
!xrandr --output DP-0 --mode 1920x1080 --rate 143.98
!xrandr --output DP-2 --mode 1920x1080 --rate 143.98
!xrandr --output DP-6 --mode 1920x1080 --rate 143.98

!xrandr --output DP-0 --mode 1920x1080 --rate 119.96
!xrandr --output DP-2 --mode 1920x1080 --rate 119.96
!xrandr --output DP-6 --mode 1920x1080 --rate 119.96

% to Manage Dual Monitor Setups
!xrandr --output DisplayPort-0 --mode 1920x1080 --rate 144.00 --output DVI-D-0 --mode 1920x1080 --rate 60.00

!xrandr --output DisplayPort-0 --primary --mode 1920x1080 --rate 144.00 --output DVI-D-0 --mode 1920x1080 --rate 60.00 --right-of DisplayPort-0

%% Set 3 monitors of the same in mirror mode
% For hnlstim2:
% The physical arrangement of the 3 monitors is as follow:
% DP-0 on the left; DP-2 (main) in the middle; DP-6 on the right;
% Don't disconnect the current DP ports! 
% Because the DP number changes every time they are reconnected!

!xrandr --output DP-0 --same-as DP-2
!xrandr --output DP-6 --same-as DP-2

%% Set 3 monitors as one in extended mode
% !xrandr --output DP-0 --same-as DP-6
% !xrandr --output DP-2 --same-as DP-6
!xrandr --output DP-2 --auto --output DP-0 --auto --left-of DP-2
!xrandr --output DP-2 --auto --output DP-6 --auto --right-of DP-2


%% Setting display port number
% first change permission to nvidia-polkit
% https://forums.developer.nvidia.com/t/cant-save-to-x-configuration-file-on-nvidia-settings/185069
!sudo chmod u+x /usr/share/screen-resolution-extra/nvidia-polkit
% Then open NVIDIA X Server Settings to make changes
% Save X configuration file to /etc/X11/xorg.conf
% Don't move the disconnect the DP ports ever again! Because the DP number
% changes every time they are disconnected.

