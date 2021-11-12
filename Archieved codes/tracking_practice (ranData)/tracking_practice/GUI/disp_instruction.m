function disp_instruction(cond_num,hText)
%this program provides instructions to the listener according to the
%condition number stored in cond_num

if cond_num == 1
    str={'In this experimental block, you will hear four sounds on each trial, one of which will be more WARBLY than the other three. Please select the one that sounds the most warbly.';'Press any key to continue.'};
elseif cond_num == 2
    str={'In this experimental block, you will hear four sounds on each trial, one of which will be more ROUGH than the other three. Please select the one that sounds the most rough.';'Press any key to continue.'};
elseif cond_num == 3
    str={'In this experimental block, you will hear four sounds on each trial, one of which will contain a GAP. Please select the one with the gap.';'Press any key to continue.'};
end
set(hText, 'String', str);
while waitforbuttonpress==0
end