!xinput  list
!xinput create-master L
!xinput  list

%% 
!xinput reattach 22 15
!xinput  list

%%
!xinput remove-master 18
!xinput  list
%%
!xinput disable 20
!xinput enable 20
%%
!xinput disable 14
!xinput enable $ID
% Disable:
!xinput set-prop 15 "Device Enabled" 0
% Enable:
!xinput set-prop 15 "Device Enabled" 1
%% set mouse speed
% make this trackball's speed twice as slow
!xinput --set-prop "Logitech USB Trackball" "Device Accel Constant Deceleration" 2
% make this trackball's speed twice as fast
!xinput --set-prop "Logitech USB Trackball" "Device Accel Constant Deceleration" 0.5

!xinput list-props 10
!xinput list-props 9
!xinput --set-prop 10 "libinput Accel Speed" 1

