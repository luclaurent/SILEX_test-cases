%%%%post processing

%file of the results
fileResult='xfem/parametricStudy/OK/2017-05-03_11-38-39_porous_cavity_400steps_ok';
%'xfem/parametricStudy/OK/2017-05-09_10-21-03_porous_cavity_90-100_100steps';
%'xfem/parametricStudy/OK/2017-05-10_16-14-55_porous_cavity_40-60_100steps';%'xfem/parametricStudy/OK/2017-05-09_10-21-03_porous_cavity_90-100_100steps';
%'xfem/parametricStudy/OK/2017-05-09_10-21-03_porous_cavity'
%'xfem/parametricStudy/OK/2017-05-10_16-14-55_porous_cavity';
%'xfem/parametricStudy/2017-05-05_18-47-45_porous_cavity';
%'xfem/parametricStudy/OK/2017-05-03_11-38-39_porous_cavity_400steps_ok';
%'xfem/parametricStudy/2017-05-02_18-41-06_porous_cavity';
%'xfem/parametricStudy/2017-05-05_12-25-22_porous_cavity_200steps_thick_ok';
%'xfem/parametricStudy/2017-05-03_11-38-39_porous_cavity_400steps_ok';
%'xfem/parametricStudy/2017-05-03_11-52-55_porous_cavity.mat';
%'classic/parametricStudy/2017-04-28_18-19-33_porous_cavity';
%'classic/parametricStudy/2017-04-28_18-03-32_porous_cavity_classic_ok_corner';
%'classic/parametricStudy/2017-04-28_18-19-33_porous_cavity';
%'classic/parametricStudy/2017-04-28_18-03-32_porous_cavity';
%'classic/parametricStudy/2017-04-28_18-19-33_porous_cavity';
%'classic/parametricStudy/2017-04-19_13-33-12_porous_cavity';
%'classic/parametricStudy/2017-04-21_19-13-19_porous_cavity';
%'classic/parametricStudy/2016-08-26_13-12-24_porous_cavity';
%'classic/parametricStudy/2016-11-25_11-56-18_porous_cavity';
fileref='Reference_noplate.mat';
addpath('matlab2tikz')
set(0,'DefaultFigureVisible','on')

%%%%%
%load results
S=load(fileResult);
R=load(fileref);

listFreqRef=R.AllFRF(1,:);

%plot all 
plotType='semilogy';
listFreq=S.varResult{1}.AllFRF(1,:);

figure;
subplot(321)
figure
for itP=1:numel(S.varResult)
    feval(plotType,listFreq,S.varResult{itP}.AllFRF(2,:));
    hold on
end
feval(plotType,listFreqRef,R.AllFRF(2,:),'k','LineWidth',2);
title('All plots')

%regroup data
funC=@(X) X.AllFRF(2,:);
tmpV=cellfun(funC,S.varResult,'UniformOutput',false);
FrFall=vertcat(tmpV{:});

%find NaN
maskNAN=isnan(FrFall);
paraNAN=S.paraValFull(maskNAN(:,end),:);
if ~isempty(paraNAN);fprintf('NAN values for parameters: \n');end
fprintf('%g %g\n',paraNAN')
%remove NaN
FrFall(maskNAN(:,end),:)=[];
% compute in dB
prefsquare=(20e-6)^2;
FrFall=10*log10(FrFall./prefsquare);



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
feval(plotType,listFreqRef,R.AllFRF(2,:),'k','LineWidth',2);
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
feval(plotType,listFreqRef,R.AllFRF(2,:),'k','LineWidth',2);
legend('CI+ (68%)+','CI- (68%)','Mean','Ref with no plate')
title('Confidence interval')
%
subplot(324)
feval(plotType,listFreq,FrFMean,'LineWidth',2,'Color','b')
hold on
errorbar(listFreq,FrFMean,FrFStd);
feval(plotType,listFreqRef,R.AllFRF(2,:),'k','LineWidth',2);
title('Error bars')
%
subplot(325)
feval(plotType,listFreq,FrFStd,'LineWidth',2,'Color','r')
hold on
feval(plotType,listFreqRef,R.AllFRF(2,:),'k','LineWidth',2);
legend('STD','Ref with no plate')

%%%% Compute mean of the FRF on a frequency window
% window
window.min=90;
window.max=100;
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

figure;
subplot(331)
for itP=1:numel(S.varResult)
    feval(plotType,listFreqW,S.varResult{itP}.AllFRF(2,ixFreqW));
    hold on
end
feval(plotType,listFreqRef(ixFreqWRef),R.AllFRF(2,ixFreqWRef),'k','LineWidth',2);
title('All plots')

%wrapped curve
subplot(432)
feval(plotType,listFreqW,FrFMin(ixFreqW),'LineWidth',2,'Color','r')
hold on
feval(plotType,listFreqW,FrFMax(ixFreqW),'LineWidth',2,'Color','r')
feval(plotType,listFreqW,FrFMean(ixFreqW),'LineWidth',2,'Color','b')
feval(plotType,listFreqRef(ixFreqWRef),R.AllFRF(2,ixFreqWRef),'k','LineWidth',2);
legend('Min','Max','Mean','Ref with no plate')
title('Min,Max,Mean')
%
plotType='plot';
%
subplot(433)
feval(plotType,listFreqW,FrFMean(ixFreqW)+FrFStd(ixFreqW),'LineWidth',2,'Color','r')
hold on
feval(plotType,listFreqW,FrFMean(ixFreqW)-FrFStd(ixFreqW),'LineWidth',2,'Color','r')
feval(plotType,listFreqW,FrFMean(ixFreqW),'LineWidth',2,'Color','b')
feval(plotType,listFreqRef(ixFreqWRef),R.AllFRF(2,ixFreqWRef),'k','LineWidth',2);
legend('CI+ (68%)+','CI- (68%)','Mean','Ref with no plate')
title('Confidence interval')
%
subplot(434)
feval(plotType,listFreqW,FrFMean(ixFreqW),'LineWidth',2,'Color','b')
hold on
errorbar(listFreqW,FrFMean(ixFreqW),FrFStd(ixFreqW));
feval(plotType,listFreqRef(ixFreqWRef),R.AllFRF(2,ixFreqWRef),'k','LineWidth',2);
title('Error bars')

subplot(435)
feval(plotType,1:nbW,FrFWMin','LineWidth',2,'Color','r')
hold on
feval(plotType,1:nbW,FrFWMax','LineWidth',2,'Color','g')
feval(plotType,1:nbW,FrFWMean','LineWidth',2,'Color','b')
%feval(plotType,1:nbW,FrFWStd','LineWidth',2,'Color','r')
feval(plotType,[1 nbW],min(R.AllFRF(2,ixFreqWRef))*ones(1,2),'k:','LineWidth',2);
feval(plotType,[1 nbW],max(R.AllFRF(2,ixFreqWRef))*ones(1,2),'k:','LineWidth',2);
feval(plotType,[1 nbW],mean(R.AllFRF(2,ixFreqWRef))*ones(1,2),'k:','LineWidth',2);
legend('Min','Max','Mean','Min no plate','Max no plate','Mean no plate')

subplot(436)
feval(plotType,listFreqW,FrFStd(ixFreqW),'LineWidth',2,'Color','r')
hold on
feval(plotType,listFreqRefW,R.AllFRF(2,ixFreqWRef),'k:','LineWidth',2);
legend('STD','STD no plate')
%
%surface
%conditionning data
XX=S.paraValFull(:,:,1);
YY=S.paraValFull(:,:,2)+1.25;
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