%% GIF of Xcorr from trial 1 to 50

for s=1:4
    
    h = figure;
    axis tight manual % this ensures that getframe() returns a consistent size
    titlename=strcat('Xcoor-theta-',conditionNames(s));
    titlename=titlename{1};
    filename =strcat(titlename,'.gif');
    

    % conditionNames={'A1&A2','A1&A1','A2&A1','A2&A2'};
    % allPerm
    % Perm=[2 3 1 4];
    % allPerm
    % conditionNames{2}


    for i=((s-1)*50+1):((s-1)*50+50)
        
        C3EEG=reRef_data(:,28,i);
        C4EEG=reRef_data(:,32,i);

        [r,lags]=xcorr(C3EEG, C4EEG);
        plot(lags./2,r);xlabel('time [ms]');ylabel('cross correlation');xlim([-250 250]);%ylim([-2*1e8 2*1e8]);
        title(titlename);

        drawnow 
        % Capture the plot as an image 
        frame = getframe(h); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        % Write to the GIF File 
        if i == ((s-1)*50+1) 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
        else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
        end 

    end
    
    clear titlename; clear filename; close all;
    
end


%% average 

for s=1:4
    subplot(4,1,s);
    titlename=strcat('Xcoor-theta-',conditionNames(s));
    titlename=titlename{1};
    
    rAll=[];

    for i=((s-1)*50+1):((s-1)*50+50)
        C3EEG=reRef_data(:,28,i);
        C4EEG=reRef_data(:,32,i);

        [r,lags]=xcorr(C3EEG, C4EEG);
        rAll=[rAll r];

    end
    rAllmean=mean(rAll,2);
    plot(lags./2,rAllmean);xlabel('time [ms]');ylabel('cross correlation');xlim([-250 250]);%ylim([-2*1e8 2*1e8]);
    title(['mean ' titlename]);
    
end
