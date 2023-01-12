figure;
tic
% cmin=0.25e-4;cmax=5e-4;
% cmin=-4e-16;cmax=4e-16;
cmin=-0.5;cmax=0.5;
for st=1:4
    pow_mat=[];
    pow_mat=squeeze(avg_state_pow(st,:,:));
    subplot(1,4,st);
    imagesc(pow_mat(1:30,:));colorbar;colormap('jet');
    % colormap(hotncold(100))
    clim([cmin cmax]);
    xlabel('chans');ylabel('freq');set(gca, 'YDir','normal');
    title(states4names{st},'color',condicolors(st,:));
end
toc
