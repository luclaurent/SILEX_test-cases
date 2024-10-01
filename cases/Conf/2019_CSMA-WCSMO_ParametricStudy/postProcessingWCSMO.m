%%%%post processing

clear all
close all

%file of the results
%fileResult='WCSMO-CSMA/paraXY800-1600.mat';
fileResult='WCSMO-CSMA/paraXY2000-1600.mat';
fileResult='WCSMO-CSMA/paraXY800-900.mat';
%fileResult='WCSMO-CSMA/paraZR1600-1600.mat';
%fileResult='WCSMO-CSMA/paraXYZ2000-1728.mat';
fileref='WCSMO-CSMA/cavity_acou3D_struc_3D_v3_nowall_results.mat';
addpath('matlab2tikz/src')
set(0,'DefaultFigureVisible','on')
%%
prefsquare=(20e-6)^2;
funFRF=@(x)10*log10(x./prefsquare);
%%%%%
%load results
S=load(fileResult);
R=load(fileref);

listFreqRef=R.AllFRF(1,:);
FRFref=funFRF(R.AllFRF(2,:));

listFreq=S.FRF(1,1,:);
listFreq=listFreq(:);
%
allFRFV=reshape(S.FRF(:,2,:),[S.nbVal^2 S.nbStep 1]);
FrFall=funFRF(allFRFV);
allFRF=cell(size(S.Xm));
for it=1:numel(S.Xm)
    allFRF{it}=funFRF(allFRFV(it,:));
end
allFRF=allFRF';
%allFRF(:,:,2:12)=[];

%plot all
plotType='plot';

figure
for itP=1:numel(allFRF)
    feval(plotType,listFreq,allFRF{itP});
    hold on
end
%feval(plotType,listFreqRef,R.AllFRF(2,:),'k','LineWidth',2);
title('All plots')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute Min, Max, Mean, Std for each frequency
FrFMin=min(FrFall);
FrFMax=max(FrFall);
FrFMean=mean(FrFall);
FrFStd=std(FrFall);
%wrapped curve
subplot(322)

feval(plotType,listFreq,FrFMin,'LineWidth',2,'Color','r')
hold on
feval(plotType,listFreq,FrFMax,'LineWidth',2,'Color','r')
feval(plotType,listFreq,FrFMean,'LineWidth',2,'Color','b')
%feval(plotType,listFreqRef,R.AllFRF(2,:),'k','LineWidth',2);
legend('Min','Max','Mean')%,'Ref with no plate')
title('Min,Max,Mean')

%
plotType='plot';
%
subplot(323)
feval(plotType,listFreq,FrFMean+FrFStd,'LineWidth',2,'Color','r')
hold on
feval(plotType,listFreq,FrFMean-FrFStd,'LineWidth',2,'Color','r')
feval(plotType,listFreq,FrFMean,'LineWidth',2,'Color','b')
%feval(plotType,listFreqRef,R.AllFRF(2,:),'k','LineWidth',2);
legend('CI+ (68%)+','CI- (68%)','Mean','Ref with no plate')
title('Confidence interval')
%
subplot(324)
feval(plotType,listFreq,FrFMean,'LineWidth',2,'Color','b')
hold on
errorbar(listFreq,FrFMean,FrFStd);
%feval(plotType,listFreqRef,R.AllFRF(2,:),'k','LineWidth',2);
title('Error bars')
%
subplot(325)
feval(plotType,listFreq,FrFStd,'LineWidth',2,'Color','r')
hold on
%feval(plotType,listFreqRef,R.AllFRF(2,:),'k','LineWidth',2);
legend('STD','Ref with no plate')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Compute mean of the FRF on a frequency window
% window

wX=linspace(S.freqMin,S.freqMax,30);
[wwX,wwY]=meshgrid(wX);

for itW=1:numel(wwX);
    window.min=wwX(itW);%40;
    window.max=wwX(itW)+wwY(itW);%100;
    if window.max < S.freqMax
        %window.min=49;
        %window.max=82;
        %frequencies in window
        ixFreqW=find(listFreq>=window.min&listFreq<=window.max);
        listFreqW=listFreq(ixFreqW);
        ixFreqWRef=find(listFreqRef>=window.min&listFreqRef<=window.max);
        listFreqRefW=listFreqRef(ixFreqWRef);
        %
        fprintf('WINDOW:    %d <f< %d\n',window.min,window.max);
        fprintf('AVAILABLE: %d <f< %d\n',min(listFreq),max(listFreq));
        fprintf('Nb steps: %i\n',numel(listFreq));
        
        %compute mean,std,min and max of the FRF in the window
        FrFWMin=min(FrFall(:,ixFreqW),[],2);
        FrFWMax=max(FrFall(:,ixFreqW),[],2);
        FrFWMean=mean(FrFall(:,ixFreqW),2);
        FrFWStd=std(FrFall(:,ixFreqW),[],2);
        
        diffMax=max(FrFWMax(:))-min(FrFWMax(:));
        diffMin=max(FrFWMin(:))-min(FrFWMin(:));
        diffMean=max(FrFWMean(:))-min(FrFWMean(:));
        fprintf('\ndiff Max:%d\n',diffMax);
        fprintf('diff Min:%d\n',diffMin);
        fprintf('diff Mean:%d\n',diffMean);
        
        
        nbW=numel(FrFWMin);
        
        %     figure;
        %     subplot(331)
        %     for itP=1:numel(allFRF)
        %         feval(plotType,listFreqW,allFRF{itP}(ixFreqW));
        %         hold on
        %     end
        %     feval(plotType,listFreqRef(ixFreqWRef),FRFref(ixFreqWRef),'k','LineWidth',2);
        %     title('All plots')
        %
        %     %wrapped curve
        %     subplot(432)
        %     feval(plotType,listFreqW,FrFMin(ixFreqW),'LineWidth',2,'Color','r')
        %     hold on
        %     feval(plotType,listFreqW,FrFMax(ixFreqW),'LineWidth',2,'Color','r')
        %     feval(plotType,listFreqW,FrFMean(ixFreqW),'LineWidth',2,'Color','b')
        %     feval(plotType,listFreqRef(ixFreqWRef),FRFref(ixFreqWRef),'k','LineWidth',2);
        %     legend('Min','Max','Mean','Ref with no plate')
        %     title('Min,Max,Mean')
        %     %
        %     plotType='plot';
        %     %
        %     subplot(433)
        %     feval(plotType,listFreqW,FrFMean(ixFreqW)+FrFStd(ixFreqW),'LineWidth',2,'Color','r')
        %     hold on
        %     feval(plotType,listFreqW,FrFMean(ixFreqW)-FrFStd(ixFreqW),'LineWidth',2,'Color','r')
        %     feval(plotType,listFreqW,FrFMean(ixFreqW),'LineWidth',2,'Color','b')
        %     feval(plotType,listFreqRef(ixFreqWRef),FRFref(ixFreqWRef),'k','LineWidth',2);
        %     legend('CI+ (68%)+','CI- (68%)','Mean','Ref with no plate')
        %     title('Confidence interval')
        %     %
        %     subplot(434)
        %     feval(plotType,listFreqW,FrFMean(ixFreqW),'LineWidth',2,'Color','b')
        %     hold on
        %     errorbar(listFreqW,FrFMean(ixFreqW),FrFStd(ixFreqW));
        %     feval(plotType,listFreqRef(ixFreqWRef),FRFref(ixFreqWRef),'k','LineWidth',2);
        %     title('Error bars')
        %
        %     subplot(435)
        %     feval(plotType,1:nbW,FrFWMin','LineWidth',2,'Color','r')
        %     hold on
        %     feval(plotType,1:nbW,FrFWMax','LineWidth',2,'Color','g')
        %     feval(plotType,1:nbW,FrFWMean','LineWidth',2,'Color','b')
        %     %feval(plotType,1:nbW,FrFWStd','LineWidth',2,'Color','r')
        %     feval(plotType,[1 nbW],min(FRFref(ixFreqWRef))*ones(1,2),'k:','LineWidth',2);
        %     feval(plotType,[1 nbW],max(FRFref(ixFreqWRef))*ones(1,2),'k:','LineWidth',2);
        %     feval(plotType,[1 nbW],mean(FRFref(ixFreqWRef))*ones(1,2),'k:','LineWidth',2);
        %     legend('Min','Max','Mean','Min no plate','Max no plate','Mean no plate')
        %
        %     subplot(436)
        %     feval(plotType,listFreqW,FrFStd(ixFreqW),'LineWidth',2,'Color','r')
        %     hold on
        %     feval(plotType,listFreqRefW,FRFref(ixFreqWRef),'k:','LineWidth',2);
        %     legend('STD','STD no plate')
        %
        
        
        
        %%plot surfaces
        
        LpMean=zeros(S.nbVal);
        LpMean(:)=mean(FrFall(:,ixFreqW),2);
        LpMean=LpMean';
        critAmp=max(LpMean(:))-min(LpMean(:));
        
        if critAmp>2
            f=figure;
            surf(S.Xm,S.Ym,LpMean)
            hlight=light;               % active light
            lighting('gouraud');         % type of rendering
            lightangle(hlight,48,70);    % direction of the light
            xlabel('$x_w$ [m]','Interpreter','latex');
            ylabel('$y_w$ [m]','Interpreter','latex');
            zlabel('$\overline{L_p}$ [dB]','Interpreter','latex');
            %ylim([0.1 1])
            %xlim([1.5,2.75])
            title(['Min ' num2str(window.min) ' Max ' num2str(window.max) ' diff ' num2str(critAmp)])
            colorbar
            %view([63 52])
            %saveas(f,['WCSMO-CSMA/Ref_obj_' num2str(window.min) '-' num2str(window.max) '.eps'],'epsc') %49-82
        end
    end
end
pause

f=figure
LpMean=zeros(S.nP);
LpMean(:)=mean(FrFallG1,2);
surf(S.Ym,S.Rm,LpMean)
hlight=light;               % active light
lighting('gouraud');         % type of rendering
lightangle(hlight,48,70);    % direction of the light
xlabel('$y_w$ [m]','Interpreter','latex');
ylabel('$R_w$ [m]','Interpreter','latex');
zlabel('$\frac{\partial\overline{L_p}}{\partial y_w}$ [dB/m]','Interpreter','latex');
colorbar
ylim([0.1 1])
xlim([1.5,2.75])
%view([63 52])
saveas(f,'WCSMO-CSMA/Ref_GX_94-102.eps','epsc')%49-82
f=figure
LpMean=zeros(S.nP);
LpMean(:)=mean(FrFallG2,2);
surf(S.Ym,S.Rm,LpMean)
hlight=light;               % active light
lighting('gouraud');         % type of rendering
lightangle(hlight,48,70);    % direction of the light
xlabel('$y_w$ [m]','Interpreter','latex');
ylabel('$R_w$ [m]','Interpreter','latex');
zlabel('$\frac{\partial\overline{L_p}}{\partial R_w}$ [dB/m]','Interpreter','latex');
colorbar
ylim([0.1 1])
xlim([1.5,2.75])
%view([63 52])
saveas(f,'WCSMO-CSMA/Ref_GY_94-102.eps','epsc')%49-82

stop

% %%%% Compute mean of the FRF on a frequency window
% % window
% window.min=90;
% window.max=100;
% %frequencies in window
% ixFreqW=find(listFreq>=window.min&listFreq<=window.max);
% listFreqW=listFreq(ixFreqW);
% ixFreqWRef=find(listFreqRef>=window.min&listFreqRef<=window.max);
% listFreqRefW=listFreqRef(ixFreqWRef);
% %
% fprintf('WINDOW:    %d <f< %d\n',window.min,window.max);
% fprintf('AVAILABLE: %d <f< %d\n',min(listFreq),max(listFreq));
% fprintf('Nb steps: %i\n',numel(listFreq));
%
% %compute mean,std,min and max of the FRF in the window
% FrFWMin=min(FrFall(:,ixFreqW),[],2);
% FrFWMax=max(FrFall(:,ixFreqW),[],2);
% FrFWMean=mean(FrFall(:,ixFreqW),2);
% FrFWStd=std(FrFall(:,ixFreqW),[],2);
%
% diffMax=max(FrFWMax(:))-min(FrFWMax(:));
% diffMin=max(FrFWMin(:))-min(FrFWMin(:));
% diffMean=max(FrFWMean(:))-min(FrFWMean(:));
% fprintf('\ndiff Max:%d\n',diffMax);
% fprintf('diff Min:%d\n',diffMin);
% fprintf('diff Mean:%d\n',diffMean);
%
%
% nbW=numel(FrFWMin);
%
% figure;
% subplot(331)
% for itP=1:numel(S.varResult)
%     feval(plotType,listFreqW,S.varResult{itP}.AllFRF(2,ixFreqW));
%     hold on
% end
% feval(plotType,listFreqRef(ixFreqWRef),R.AllFRF(2,ixFreqWRef),'k','LineWidth',2);
% title('All plots')
%
% %wrapped curve
% subplot(432)
% feval(plotType,listFreqW,FrFMin(ixFreqW),'LineWidth',2,'Color','r')
% hold on
% feval(plotType,listFreqW,FrFMax(ixFreqW),'LineWidth',2,'Color','r')
% feval(plotType,listFreqW,FrFMean(ixFreqW),'LineWidth',2,'Color','b')
% feval(plotType,listFreqRef(ixFreqWRef),R.AllFRF(2,ixFreqWRef),'k','LineWidth',2);
% legend('Min','Max','Mean','Ref with no plate')
% title('Min,Max,Mean')
% %
% plotType='plot';
% %
% subplot(433)
% feval(plotType,listFreqW,FrFMean(ixFreqW)+FrFStd(ixFreqW),'LineWidth',2,'Color','r')
% hold on
% feval(plotType,listFreqW,FrFMean(ixFreqW)-FrFStd(ixFreqW),'LineWidth',2,'Color','r')
% feval(plotType,listFreqW,FrFMean(ixFreqW),'LineWidth',2,'Color','b')
% feval(plotType,listFreqRef(ixFreqWRef),R.AllFRF(2,ixFreqWRef),'k','LineWidth',2);
% legend('CI+ (68%)+','CI- (68%)','Mean','Ref with no plate')
% title('Confidence interval')
% %
% subplot(434)
% feval(plotType,listFreqW,FrFMean(ixFreqW),'LineWidth',2,'Color','b')
% hold on
% errorbar(listFreqW,FrFMean(ixFreqW),FrFStd(ixFreqW));
% feval(plotType,listFreqRef(ixFreqWRef),R.AllFRF(2,ixFreqWRef),'k','LineWidth',2);
% title('Error bars')
%
% subplot(435)
% feval(plotType,1:nbW,FrFWMin','LineWidth',2,'Color','r')
% hold on
% feval(plotType,1:nbW,FrFWMax','LineWidth',2,'Color','g')
% feval(plotType,1:nbW,FrFWMean','LineWidth',2,'Color','b')
% %feval(plotType,1:nbW,FrFWStd','LineWidth',2,'Color','r')
% feval(plotType,[1 nbW],min(R.AllFRF(2,ixFreqWRef))*ones(1,2),'k:','LineWidth',2);
% feval(plotType,[1 nbW],max(R.AllFRF(2,ixFreqWRef))*ones(1,2),'k:','LineWidth',2);
% feval(plotType,[1 nbW],mean(R.AllFRF(2,ixFreqWRef))*ones(1,2),'k:','LineWidth',2);
% legend('Min','Max','Mean','Min no plate','Max no plate','Mean no plate')
%
% subplot(436)
% feval(plotType,listFreqW,FrFStd(ixFreqW),'LineWidth',2,'Color','r')
% hold on
% feval(plotType,listFreqRefW,R.AllFRF(2,ixFreqWRef),'k:','LineWidth',2);
% legend('STD','STD no plate')
% %
% %surface
% %conditionning data
XX=S.paraValFull(:,:,1);
YY=S.paraValFull(:,:,2);
%display data
dispZ='FrFWMean';
nbLvl=40;
ZZ=reshape(eval(dispZ),size(XX));
subplot(437)
%surf(XX,YY,mean(R.AllFRF(2,ixFreqWRef))*ones(size(XX)))
%alpha(.4)
%hold on
surf(XX,YY,ZZ)
zlabel(dispZ)
subplot(4,3,10)
zmin = floor(min(ZZ(:)));
zmax = ceil(max(ZZ(:)));
zinc = (zmax - zmin) / nbLvl;
zlevs = zmin:zinc:zmax;
contour(XX,YY,ZZ,zlevs)


dispZ='FrFWMin';
ZZ=reshape(eval(dispZ),size(XX));
subplot(4,3,8)
surf(XX,YY,ZZ)
zlabel(dispZ)
subplot(4,3,11)
zmin = floor(min(ZZ(:)));
zmax = ceil(max(ZZ(:)));
zinc = (zmax - zmin) / nbLvl;
zlevs = zmin:zinc:zmax;
contour(XX,YY,ZZ,zlevs)


dispZ='FrFWMax';
ZZ=reshape(eval(dispZ),size(XX));
subplot(4,3,9)
surf(XX,YY,ZZ)
zlabel(dispZ)
subplot(4,3,12)
zmin = floor(min(ZZ(:)));
zmax = ceil(max(ZZ(:)));
zinc = (zmax - zmin) / nbLvl;
zlevs = zmin:zinc:zmax;
contour(XX,YY,ZZ,zlevs)





%export
f=figure;
dispZ='FrFWMean';
ZZ=reshape(eval(dispZ),size(XX));
s=surf(XX,YY,ZZ)
hlight=light;               % active light
lighting('gouraud');         % type of rendering
lightangle(hlight,48,70);    % direction of the light
%set(s,'EdgeColor','w')
set(s,'LineWidth',1)
colorbar
view(-63,59)
%hold on
%contour(XX,YY,ZZ,zlevs)
%hold off
xlabel('$\theta$','Interpreter','latex')
ylabel('$d$','Interpreter','latex')
zlabel('$F_{obj}$','Interpreter','latex')

saveas(f,'export/Ref_objm.pdf')
saveas(f,'export/Ref_objm.eps','epsc')
%matlab2tikz('standalone',true,'export/Ref_obj.tex')

% %contourf(XX,YY,ZZ,zlevs,'ShowText','on')
f=figure;
nbLvl=10;
zmin = floor(min(ZZ(:)));
zmax = ceil(max(ZZ(:)));
zinc = (zmax - zmin) / nbLvl;
zlevs = zmin:zinc:zmax;
[C,h]=contourf(XX,YY,ZZ,10,'ShowText','on');
clabel(C,h,'BackgroundColor',[1 1 .6],'Edgecolor',[.7 .7 .7]);
set(h,'LineWidth',2);
xlabel('$\theta$','Interpreter','latex')
ylabel('$d$','Interpreter','latex')
colorbar

saveas(f,'export/Ref_obj_contourm.pdf')
saveas(f,'export/Ref_obj_contourm.eps','epsc')
%matlab2tikz('standalone',true,'export/Ref_obj_contour.tex')