%%%%post processing

%file of the results
fileResult='2016-11-25_15-39-51_runDOE.mat';
%fileref='Reference_noplate.mat';

%%%%%
%load results
S=load(fileResult);

%%%% Compute mean of the FRF on a frequency window
% window
window.min=43;
window.max=62;

%
nP=30;
[XX,YY]=meshgrid(linspace(S.angle.min,S.angle.max,nP),linspace(S.position.min,S.position.max,nP));
gridRef(:,:,1)=XX;
gridRef(:,:,2)=YY;

itSILEX=1;
samplePts=S.SILEX(itSILEX).paraValFull;
for itS=1:size(samplePts,1);
    resp(itS)=mean(S.SILEX(1).varResult{itSILEX}{itS}.AllFRF(2,:));
end


addpath('grenat')

Meta=GRENAT('KRG',samplePts,resp);

Meta.defineRef(gridRef);
Meta.show;




pause

%plot all 
plotType='semilogy';
listFreq=S.varResult{1}.AllFRF(1,:);

figure;
subplot(321)
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
legend('Min','Max','Mean','Ref with no plate')
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
window.min=43;
window.max=62;
%frequencies in window
ixFreqW=find(listFreq>=window.min&listFreq<=window.max);
listFreqW=listFreq(ixFreqW);
ixFreqWRef=find(listFreqRef>=window.min&listFreqRef<=window.max);
listFreqRefW=listFreqRef(ixFreqWRef);

%compute mean,std,min and max of the FRF in the window
FrFWMin=min(FrFall(:,ixFreqW),[],2);
FrFWMax=max(FrFall(:,ixFreqW),[],2);
FrFWMean=mean(FrFall(:,ixFreqW),2);
FrFWStd=std(FrFall(:,ixFreqW),[],2);
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
YY=S.paraValFull(:,:,2);
%display data
dispZ='FrFWMean';
nbLvl=40;
ZZ=reshape(eval(dispZ),size(XX));
subplot(437)
surf(XX,YY,mean(R.AllFRF(2,ixFreqWRef))*ones(size(XX)))
alpha(.4)
hold on
surf(XX,YY,ZZ)
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
subplot(4,3,12)
zmin = floor(min(ZZ(:)));
zmax = ceil(max(ZZ(:)));
zinc = (zmax - zmin) / nbLvl;
zlevs = zmin:zinc:zmax;
contour(XX,YY,ZZ,zlevs)


%matlab2tikz('standalone',true,'export/Ref_obj.tex')

% %contourf(XX,YY,ZZ,zlevs,'ShowText','on')
% [C,h]=contourf(XX,YY,ZZ,zlevs,'ShowText','on')
% clabel(C,h,'BackgroundColor',[1 1 .6],'Edgecolor',[.7 .7 .7]);
% set(h,'LineWidth',2);
% xlabel('\theta')
% ylabel('d')
% matlab2tikz('standalone',true,'export/Ref_obj_contour.tex')