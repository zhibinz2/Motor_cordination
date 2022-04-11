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

% to Manage Dual Monitor Setups
!xrandr --output DisplayPort-0 --mode 1920x1080 --rate 144.00 --output DVI-D-0 --mode 1920x1080 --rate 60.00

!xrandr --output DisplayPort-0 --primary --mode 1920x1080 --rate 144.00 --output DVI-D-0 --mode 1920x1080 --rate 60.00 --right-of DisplayPort-0

%% Set 3 monitors in mirror mode
!xrandr --output DP-0 --same-as DP-2
!xrandr --output DP-6 --same-as DP-2