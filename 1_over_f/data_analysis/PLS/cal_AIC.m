%% Try AIC
% open /home/zhibin/Documents/sedation-restingstate/20200522CSP/classifyRS_better_version.m

% Calculate AIC
AIC(i)=2*log(ERRORRATE)+2*(ncomps(i));

ssEr

R2

%% Try mynpls_pred ver 1
% clear plsmodel;
% cmin=-8e-3;cmax=8e-3;
% cmin=-8e-3;cmax=8e-3;
AIC=[];c=[1 4];R2s=[];
for Fac=1:10;
    for s=1:2 % one of the four states
        R2=[];reg=[];ypred=[];ssEr=[];
        [R2,reg,ypred,ssEr] = mynpls_pred(pow5forpls3(Inds4(:,c(s)),:),H_all_LR(Inds4(:,c(s))),Fac);
%         plsmodel(c).weights=reshape(reg{Fac},5,32);
%         plsmodel(c).R2 = R2;
%         plsmodel(c).ypred = ypred;
%         subplot(2,2,c);
%         imagesc(plsmodel(c).weights);colorbar; 
%         yticks([1:5]);yticklabels({'delta','theta','alpha','beta','gamma'});
%         set(gca, 'YDir','normal');
%         xticks([1:32]);xticklabels([labels]);xtickangle(90);
%         colormap('jet'); clim([cmin cmax]);% caxis([-0.1 0.1]);
%         title([states4names{c} ': PLS model (R2= ' num2str(round(R2,1)) '  Fac= ' num2str(Fac) ') '], ...
%             'Color',condicolors(c,:));
%         grid on;
        R2s(s,Fac)=R2;
        % Cal AIC
        AIC(s,Fac)=48*log(ssEr)+2*(Fac);
    end
end

figure;
subplot(2,2,1);
plot(1:10,AIC(1,:),'r.');xlabel('nFac');ylabel('AIC');title('AIC(Fac)=48*log(ssEr)+2*(Fac)')
subtitle('uncouple')
subplot(2,2,2)
plot(1:10,AIC(2,:),'r.');xlabel('nFac');ylabel('AIC');title('AIC(Fac)=48*log(ssEr)+2*(Fac)')
subtitle('mutual')
subplot(2,2,3);
plot(1:10,R2s(1,:),'r.');xlabel('nFac');ylabel('R2');
subtitle('uncouple')
subplot(2,2,4)
plot(1:10,R2s(2,:),'r.');xlabel('nFac');ylabel('R2');
subtitle('mutual')

%% Try mynpls_pred ver 2 (commented out cross validation )
% clear plsmodel;
% cmin=-8e-3;cmax=8e-3;
% cmin=-8e-3;cmax=8e-3;
c=[1 4];
AIC=nan(2,10);R2s=nan(2,10);ssERs=nan(2,10);
    for s=1:2 % one of the four states
        reg=[];X=[];Y=[];
        X=pow5forpls3(Inds4(:,c(s)),:);Y=H_all_LR(Inds4(:,c(s)));
        [reg,ypred_fit,X_MC,Y_MC] = mynpls_pred(X,Y);
        ssOr=sum(sum((Y-mean(Y)).^2))
        for Fac=1:10;
            R2=[];ssEr=[];
            ssEr=sum(sum((Y-ypred_fit(:,Fac)).^2))% ssEr=sum(sum((Y-ypred).^2))
            ssERs(s,Fac)=ssEr;
            R2=100*(1-ssEr/ssOr)
            R2s(s,Fac)=R2;
            % Cal AIC
            AIC(s,Fac)=log(ssEr)+2*(Fac);
    end
end
figure;
subplot(3,2,1);
plot(1:10,AIC(1,:),'r.');xlabel('nFac');ylabel('AIC');title('AIC(s,Fac)=log(ssEr)+2*(Fac);')
subtitle('uncouple')
subplot(3,2,2)
plot(1:10,AIC(2,:),'r.');xlabel('nFac');ylabel('AIC');title('AIC(s,Fac)=log(ssEr)+2*(Fac);')
subtitle('mutual')
subplot(3,2,3);
plot(1:10,R2s(1,:),'r.');xlabel('nFac');ylabel('R2');
subtitle('uncouple')
subplot(3,2,4)
plot(1:10,R2s(2,:),'r.');xlabel('nFac');ylabel('R2');
subtitle('mutual')
set(gcf,'color','w'); % set background white for copying in ubuntu
subplot(3,2,5);
plot(1:10,ssERs(1,:),'r.');xlabel('nFac');ylabel('ssEr');% title('AIC(s,Fac)=log(ssEr)+2*(Fac);')
subtitle('uncouple')
subplot(3,2,6)
plot(1:10,ssERs(2,:),'r.');xlabel('nFac');ylabel('ssEr');% title('AIC(s,Fac)=log(ssEr)+2*(Fac);')
subtitle('mutual')
%% Mean center everything (This is correct)
c=[1 4];regs=[];
AIC=nan(2,10);R2s=nan(2,10);ssERs=nan(2,10);
    for s=1:2 % one of the four states
        reg=[];X=[];Y=[];
        X=pow5forpls3(Inds4(:,c(s)),:);Y=H_all_LR(Inds4(:,c(s)));
        [reg,ypred_fit,X_MC,Y_MC] = mynpls_pred(X,Y);
        regs{s}=reg;
        ssOr=sum(sum((Y_MC-mean(Y_MC)).^2))
        for Fac=1:10;
            R2=[];ssEr=[];
            ssEr=sum(sum((Y_MC-ypred_fit(:,Fac)).^2))
            ssERs(s,Fac)=ssEr;
            R2=100*(1-ssEr/ssOr)
            R2s(s,Fac)=R2;
            % Cal AIC
            AIC(s,Fac)=log(ssEr)+2*(Fac);
    end
end
figure;
subplot(3,2,1);
plot(1:10,AIC(1,:),'r.');xlabel('nFac');ylabel('AIC');title('AIC(s,Fac)=log(ssEr)+2*(Fac);')
subtitle('uncouple')
subplot(3,2,2)
plot(1:10,AIC(2,:),'r.');xlabel('nFac');ylabel('AIC');title('AIC(s,Fac)=log(ssEr)+2*(Fac);')
subtitle('mutual')
subplot(3,2,3);
plot(1:10,R2s(1,:),'r.');xlabel('nFac');ylabel('R2');
subtitle('uncouple')
subplot(3,2,4)
plot(1:10,R2s(2,:),'r.');xlabel('nFac');ylabel('R2');
subtitle('mutual')
set(gcf,'color','w'); % set background white for copying in ubuntu
subplot(3,2,5);
plot(1:10,ssERs(1,:),'r.');xlabel('nFac');ylabel('ssEr');% title('AIC(s,Fac)=log(ssEr)+2*(Fac);')
subtitle('uncouple')
subplot(3,2,6)
plot(1:10,ssERs(2,:),'r.');xlabel('nFac');ylabel('ssEr');% title('AIC(s,Fac)=log(ssEr)+2*(Fac);')
subtitle('mutual')

%

regs{1}{1}

clear plsmodel
% reshape regs{1}{1} and {2}{1}
plsmodel(1).weights=reshape(regs{1}{Fac},5,32); % uncouple
plsmodel(2).weights=reshape(regs{2}{Fac},5,32); % mutual

% topoplot for uncouple and mutual (same as in PLOT 16)
addpath /home/zhibin/Documents/GitHub/matlab-archive/hnlcode/common/gen_code/color
hnc = hotncold(100);
band5names={'delta','theta','alpha','beta','gamma'};
states2names={'uncouple','mutual'};
states4names;
figure; % canvas(0.9,0.5)
cmin=-0.01;cmax=0.01;
c=[1 4];
for s=1:2
    for b=1:5
        subplot(2,5,(s-1)*5+b)
        topoplot(plsmodel(s).weights(b,:),channels,'nosedir','+X');colorbar; 
        colormap(hnc); % colormap('jet'); 
        clim([cmin cmax]);
        if b==3;title({['PLS model: sum-EEG(-500ms) -> H-int ^{* PLOT 16}'], ...
                [states4names{c(s)} '(R^2= ' num2str(round(R2s(s),1)) ...
            '  Fac= ' num2str(Fac) ') ']},'Color',condicolors(c(s),:));end
        subtitle(band5names(b));
    end
end
set(gcf,'color','w'); % set background white for copying in ubuntu
