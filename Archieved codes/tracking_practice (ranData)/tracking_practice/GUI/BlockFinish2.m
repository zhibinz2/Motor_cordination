% BlockFinish2.m
%taken from block finish, but does not revers the axis

data.labdata.avg = levitt('threshold',1);
data.labdata.sd = levitt('deviation',1);
set(hPlot, 'XData', 1:data.param.ntrials );
set(hPlot,'YData', data.labdata.slev);
%set(hAxes, 'YDir', 'reverse');
sresult = sprintf('Mean is: %5.2f,  Stdev is: %5.2f', data.labdata.avg, data.labdata.sd);
save(['..\data\' data.param.rootfile '\' data.param.rootfile int2str(data.param.itfile)],'data');

set(hTitle,'String',sresult);
set(hPlot,'Visible','on');
set(hAxes,'Visible','on');
set(hText,'Visible','off');
while waitforbuttonpress == 0
end
set(hAxes,'Visible','off');
set(hPlot,'Visible','off');
set(hText,'Visible','on');
