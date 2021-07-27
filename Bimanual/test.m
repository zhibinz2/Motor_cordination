open TouchInput

TouchDevID=GetTouchDeviceIndices()
GetTouchDeviceInfo(TouchDevID)
Create (TouchQueueCreate()) 

TouchQueueStart()
TouchQueueStop()
TouchEventFlush()
% to start and stop data collection, or to discard all currently pending touch input.

TouchEventAvail() %to find the number of collected touch events

TouchEventGet() %to retrieve events. These events provide all the basic
%    information of interest: The time when the event was received, the (x,y)
%    touch coordinates in different useful coordinate systems, if it was a
%    start of a new touch, the end of an old touch, or some state change on an
%    existing touch, e.g., position, pressure, shape etc. Additional device
%    and operating system specific info (non-portable!) can be retrieved via
%    the helper function GetTouchValuators().


[event, nremaining] = TouchEventGet(deviceIndex, windowHandle [, maxWaitTimeSecs=0])

screenId=0;
MultiTouchMinimalDemo(TouchDevID,screenId);


%% https://wiki.archlinux.org/title/Multi-pointer_X

% zhibin@hnlstim1-ss-uci-edu:~$ xinput list
% ⎡ Virtual core pointer                    	id=2	[master pointer  (3)]
% ⎜   ↳ Virtual core XTEST pointer              	id=4	[slave  pointer  (2)]
% ⎜   ↳ Logitech USB Optical Mouse              	id=13	[slave  pointer  (2)]
% ⎜   ↳ Wacom Intuos Pro L Pad pad              	id=9	[slave  pointer  (2)]
% ⎜   ↳ Wacom Intuos Pro L Pen stylus           	id=10	[slave  pointer  (2)]
% ⎜   ↳ Wacom Intuos Pro L Pen eraser           	id=11	[slave  pointer  (2)]
% ⎜   ↳ Wacom Intuos Pro L Pen cursor           	id=14	[slave  pointer  (2)]
% ⎜   ↳ Wacom Intuos Pro L Finger touch         	id=15	[slave  pointer  (2)]
% ⎜   ↳ USB Optical Mouse                       	id=16	[slave  pointer  (2)]
% ⎣ Virtual core keyboard                   	id=3	[master keyboard (2)]
%     ↳ Virtual core XTEST keyboard             	id=5	[slave  keyboard (3)]
%     ↳ Power Button                            	id=6	[slave  keyboard (3)]
%     ↳ Power Button                            	id=7	[slave  keyboard (3)]
%     ↳ Sleep Button                            	id=8	[slave  keyboard (3)]
%     ↳ Dell Dell USB Keyboard                  	id=12	[slave  keyboard (3)]

%%
zhibin@hnlstim1-ss-uci-edu:~$ xinput list
% ⎡ Virtual core pointer                    	id=2	[master pointer  (3)]
% ⎜   ↳ Virtual core XTEST pointer              	id=4	[slave  pointer  (2)]
% ⎜   ↳ Logitech USB Optical Mouse              	id=13	[slave  pointer  (2)]
% ⎜   ↳ Wacom Intuos Pro L Finger touch         	id=15	[slave  pointer  (2)]
% ⎜   ↳ Wacom Intuos Pro L Pad pad              	id=9	[slave  pointer  (2)]
% ⎜   ↳ Wacom Intuos Pro L Pen stylus           	id=10	[slave  pointer  (2)]
% ⎜   ↳ Wacom Intuos Pro L Pen eraser           	id=11	[slave  pointer  (2)]
% ⎜   ↳ Wacom Intuos Pro L Pen cursor           	id=14	[slave  pointer  (2)]
% ⎣ Virtual core keyboard                   	id=3	[master keyboard (2)]
%     ↳ Virtual core XTEST keyboard             	id=5	[slave  keyboard (3)]
%     ↳ Power Button                            	id=6	[slave  keyboard (3)]
%     ↳ Power Button                            	id=7	[slave  keyboard (3)]
%     ↳ Sleep Button                            	id=8	[slave  keyboard (3)]
%     ↳ Dell Dell USB Keyboard                  	id=12	[slave  keyboard (3)]
% ⎡ Auxiliary pointer                       	id=17	[master pointer  (18)]
% ⎜   ↳ USB Optical Mouse                       	id=16	[slave  pointer  (17)]
% ⎜   ↳ Auxiliary XTEST pointer                 	id=19	[slave  pointer  (17)]
% ⎣ Auxiliary keyboard                      	id=18	[master keyboard (17)]
%     ↳ Auxiliary XTEST keyboard                	id=20	[slave  keyboard (18)]


%%
clear
DuoMiceTrace
plot(theX,theY);


%%   **********  https://psychtoolbox.discourse.group/t/wacom-tablet-valuators/3393
while drawing 
        [x,y,touch,~,val{seg},valinfo] = GetMouse(wdwprop.window); 
        t=GetSecs;

        if ~touch(1) %stopped drawing
            drawing=false;
        
        else % still drawing, continue line
            
            if x ~= Xpos || y ~= Ypos %store if mouse moved
                time{seg}= [time{seg} t-start];
                points{seg} = [points{seg} ; x y]; 
                [numPoints, ~]=size(points{seg});

            % draw the most recent line segment
                if feedback
                    Screen(wdwprop.window,'DrawLine',128,points{seg}(numPoints-1,1),points{seg}(numPoints-1,2),points{seg}(numPoints,1),points{seg}(numPoints,2));
                    % add new segment without clearing, don't wait for refresh reduces lag
                    Screen('Flip', wdwprop.window, 0, 1, 2);
                end
                Xpos=x; Ypos=y;
            end
        end                
end
    
%%
[x,y,touch,~,val,valinfo] = GetMouse(wdwprop.window,7); 
%%
name = 
    {'Virtual core pointer'          } 
    {'Virtual core XTEST pointer'    }
    {'Wacom Cintiq 22HD Pen stylus'  }
    {'Wacom Cintiq 22HD Pad pad'     }
    {'DELL09BD:00 0488:120A Touchpad'}
    {'PS/2 Generic Mouse'            }
    {'Wacom Cintiq 22HD Pen eraser'  }

info{3}=
  struct with fields: 

     usagePageValue: 1
         usageValue: 2
          usageName: 'slave pointer'
              index: 7
          transport: 'enabled'
           vendorID: []
          productID: []
            version: []
       manufacturer: []
            product: 'Wacom Cintiq 22HD Pen stylus'
       serialNumber: []
         locationID: 2
        interfaceID: 9
      totalElements: 262
           features: 8
             inputs: 262
            outputs: 0
        collections: 0
               axes: 6
            buttons: 256
               hats: 0
            sliders: 0
              dials: 0
             wheels: 0
    touchDeviceType: -1
     maxTouchpoints: -1

