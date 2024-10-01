%%%%post processing 
%% check mesh convergence

%parameters: 180 2.7

%file of the results
fileResult{1}='2017-04-19_13-38-19_porous_cavity.mat';
critMesh(1)=80;
fileResult{2}='2017-04-19_13-33-12_porous_cavity.mat';
critMesh(2)=60;
fileResult{3}='2017-04-19_13-30-37_porous_cavity.mat';
critMesh(3)=40;
fileResult{4}='2017-04-19_13-29-22_porous_cavity.mat';
critMesh(4)=20;
fileref='Reference_noplate.mat';

%%%%%
%load results
S={};
for itL=1:numel(fileResult)
    S{itL}=load(fileResult{itL});
end
R=load(fileref);

listFreqRef=R.AllFRF(1,:);

%plot all 
plotType='semilogy';
listFreq={};
FRF={};
for itL=1:numel(S)
    listFreq{itL}=S{itL}.varResult{1}.AllFRF(1,:);
    FRF{itL}=S{itL}.varResult{1}.AllFRF(2,:);
end

figure;
for itL=1:numel(S)
    feval(plotType,listFreq{itL},FRF{itL});
    hold on
end
pause



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
window.min=50;
window.max=60;
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
YY=S.paraValFull(:,:,2)+1.25;
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


%export
f=figure
dispZ='FrFWMean';
ZZ=reshape(eval(dispZ),size(XX));
surf(XX,YY,ZZ)
xlabel('$\theta$','Interpreter','latex')
ylabel('$d$','Interpreter','latex')
zlabel('$F_{obj}$','Interpreter','latex')
saveas(f,'export/Ref_objm.pdf')
saveas(f,'export/Ref_objm.eps')
matlab2tikz('standalone',true,'export/Ref_obj.tex')

% %contourf(XX,YY,ZZ,zlevs,'ShowText','on')
f=figure
[C,h]=contourf(XX,YY,ZZ,10,'ShowText','on')
clabel(C,h,'BackgroundColor',[1 1 .6],'Edgecolor',[.7 .7 .7]);
set(h,'LineWidth',2);
xlabel('$\theta$','Interpreter','latex')
ylabel('$d$','Interpreter','latex')
saveas(f,'export/Ref_obj_contourm.pdf')
saveas(f,'export/Ref_obj_contourm.eps')
matlab2tikz('standalone',true,'export/Ref_obj_contour.tex')