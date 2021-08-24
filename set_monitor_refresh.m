% to Find Info About Your Monitors
!xrandr -q

% to Set Your Monitor Resolution
!xrandr --output DVI-D-0 --mode 1920x1080

% to Set Your Monitor Refresh Rate
!xrandr --output DVI-D-0 --mode 1920x1080 --rate 60.00

% to Manage Dual Monitor Setups
!xrandr --output DisplayPort-0 --mode 1920x1080 --rate 144.00 --output DVI-D-0 --mode 1920x1080 --rate 60.00

!xrandr --output DisplayPort-0 --primary --mode 1920x1080 --rate 144.00 --output DVI-D-0 --mode 1920x1080 --rate 60.00 --right-of DisplayPort-0