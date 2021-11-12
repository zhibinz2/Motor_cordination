function y=fdbk(cor,hText)
%this program provides feedback for command window
%  and returns correct (1) or incorrect (0)
%tempStr = '            real answer       user answer';
%tempStr1 =['                     ',int2str(ran),'       ',int2str(answer)];
if cor==1
   %   tempStr2='                       Right ';
   str=' Right ';
   %   str=str2mat(tempStr, tempStr1, ' ', tempStr2);
   y=1;
else
   str=' Wrong ';
   %   tempStr2='                       wrong ';
   %   str=str2mat(tempStr, tempStr1, ' ', tempStr2);
   y=0;
end;
set(hText, 'String', str);
drawnow;
pause(1);
set(hText, 'String', ' ');
pause(0.5);