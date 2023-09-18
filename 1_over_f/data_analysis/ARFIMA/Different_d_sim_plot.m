%%  color scheme
red   = [1 0 0];
pink  = [1 0.65 0.75];
blue  = [0 0 1];
green = [0 1 0];
darkgreen = [0 0.5 0];
grey  = [0.5 0.5 0.5];
yellow  = [1 1 0];
deepyellow  = [1 0.8 0.2];
megenta = [1 0 1];% fill([0 1 1 0],[0 0 1 1],megenta)
cyan = [0 1 1]; % fill([0 1 1 0],[0 0 1 1],cc)
purple = [0.6 0.1 0.9];
ARFIMAcolors=[deepyellow;darkgreen;pink];

%% PLOT-3 Different d
figure('units','normalized','outerposition',[0 0 1 1]);
clf
N=1000;stdx=20;
ds=[0.2 0.5 1];
for i=1:length(ds)
    [Z] = ARFIMA_SIM(N,[],[],ds(i),stdx);
    y=Z+750;
    subplot(length(ds),4,(i-1)*4+1); 
        plot(y,'b');
%         title(['d = ' num2str(ds(i))]);
        ax = gca;
        ax.FontSize = 13; 
        xlabel('Taps');ylabel('Tapping interval (ms)')
        ylim([680 840]);
    subplot(length(ds),4,(i-1)*4+2); 
        autocorr(y,N-1);
%         title(['d = ' num2str(ds(i))]);
        title([]);
        ylabel('Autocorrelation (\rho)')
        ax = gca;
        ax.FontSize = 13; 
        ylim([-0.5 1]);
    subplot(length(ds),4,(i-1)*4+3); 
        Fs=1.3; 
        [freqs,fcoef,beta,xx,yy,FitValues] = oneoverf(y,Fs);
%         title(['d = ' num2str(ds(i)) '   (Beta = d x 2)']);
        ax = gca;
        ax.FontSize = 13; 
        ylim([-6 8]);
    subplot(length(ds),4,(i-1)*4+4); 
        [D,Alpha1,n,F_n,FitValues]=DFA_main(y);
%         title(['d = ' num2str(ds(i)) '   (H = d+0.5)']);
        ax = gca;
        ax.FontSize = 13; 
        ylim([0 3]);
end
set(gcf,'color','w'); % set background white for copying in ubuntu
sgtitle(['Different d with std = ' num2str(stdx)]);

%%
Ns=[100 600 1000]; ds=[0.1:0.1:1]; repeat=10; stdx=20;
figure('units','normalized','outerposition',[0 0 0.25 0.45]);

clf

for i=1:3
    [Betass,BetassErr,Hss,HssErr] = IterARFIMA(Ns(i),ds,repeat,stdx);

    subplot(1,2,1);
        hold on;
        errorbar(ds,mean(Betass,2),BetassErr,'.','color',ARFIMAcolors(i,:));

        hold off;
        if i==3
            B=2.*ds; hold on; plot(ds,B,'color',purple);hold off;
            legend({'n = 100','n = 600','n = 1000','Beta = 2 x d'},'Location','southeast');
            xlim([0 1.1]);ylim([0 2.5])
            xlabel('d');ylabel('Beta');
            title(['Spectral Analysis']);
            ax = gca;
            ax.FontSize = 17; 
        end
    subplot(1,2,2)
        hold on;
        errorbar(ds,mean(Hss,2),HssErr,'.','color',ARFIMAcolors(i,:))
        
        hold off;
        if i==3
            H=ds+0.5; hold on; plot(ds,H,'color',purple);hold off;
            legend({'n = 100','n = 600','n = 1000','H = d + 0.5'},'Location','southeast');
            xlim([0 1.1]);ylim([0.5 1.6])
            xlabel('d');ylabel('H');
            title(['DFA']);
            ax = gca;
            ax.FontSize = 17; 
        end
end

set(gcf,'color','w'); % set background white for copying in ubuntu