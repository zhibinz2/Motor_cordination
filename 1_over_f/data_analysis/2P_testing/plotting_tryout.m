clc,clear
clf
%%
figure
subplot 121
plot(0,0,'r.')
subplot 122
plot(1,1,'b.')
h1 = annotation('textbox',[0.25 0.5 0.1 0.05],...
    'string','PLOT1');
h2 = annotation('textbox',[0.65 0.5 0.1 0.05],...
    'string','PLOT2');
set([h1 h2], 'fitboxtotext','on',...
    'edgecolor','none')


%%
figure
% dim is [x position, y position of starting point, width, height]
h0=annotation('textbox',[0/5+2/20 0/5+19/20 1/5/2 1/20],'string','A','color',red)
h1=annotation('textbox',[1/5+2/20 0/5+19/20 1/5/2 1/20],'string','B','color',red)
h2=annotation('textbox',[2/5+2/20 0/5+19/20 1/5/2 1/20],'string','C','color',red)
h3=annotation('textbox',[3/5+2/20 0/5+19/20 1/5/2 1/20],'string','D','color',red)
h4=annotation('textbox',[4/5+2/20 0/5+19/20 1/5/2 1/20],'string','E','color',red)
v0=annotation('textbox',[0/5+1/20 0/5+1/20 1/5/2 1/20],'string','H','color',red)
v1=annotation('textbox',[0/5+1/20 1/5+1/20 1/5/2 1/20],'string','I','color',red)
v2=annotation('textbox',[0/5+1/20 2/5+1/20 1/5/2 1/20],'string','J','color',red)
v3=annotation('textbox',[0/5+1/20 3/5+1/20 1/5/2 1/20],'string','K','color',red)
set([h0 h1 h2 h3 h4 v0 v1 v2 v3], 'fitboxtotext','on',...
    'edgecolor','none')


%%
https://www.mathworks.com/help/matlab/ref/tiledlayout.html

https://www.mathworks.com/help/matlab/ref/colorbar.html

https://www.mathworks.com/help/matlab/ref/matlab.graphics.illustration.colorbar-properties.html

https://www.mathworks.com/matlabcentral/answers/60376-how-to-make-an-inset-of-matlab-figure-inside-the-figure