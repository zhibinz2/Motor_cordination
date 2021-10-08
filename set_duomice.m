!xinput  list
!xinput create-master L
!xinput  list

%% dell
!xinput reattach 19 20
!xinput reattach 9 20
!xinput  list

%% BenQ

!xinput reattach 9 16
!xinput reattach 13 2
!xinput reattach 15 9
!xinput  list
!xinput remove-master 20

!xinput  list
!xinput create-master L
!xinput  list
!xinput reattach 12 14
%% hnlb
!xinput reattach 14 16
!xinput reattach 29 16
!xinput reattach 12 14
!xinput reattach 9 16
!xinput reattach 12 9
!xinput  list

%% set mouse speed
% make this trackball's speed twice as slow
!xinput --set-prop "Logitech USB Trackball" "Device Accel Constant Deceleration" 2


% make this trackball's speed twice as fast
!xinput --set-prop "Logitech USB Trackball" "Device Accel Constant Deceleration" 0.5


!xinput list-props 10
!xinput list-props 9
!xinput --set-prop 10 "libinput Accel Speed" 1

