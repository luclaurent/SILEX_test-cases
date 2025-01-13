%% execute SILEX for building kriging surrogate model
clear all
close all
%add path
addpath('grenat')

%load folder structure
initDirGRENAT;
%
nsmin=2;
nsmax=10;
nbRep=1;
nscomb=nsmin:nsmax;
%compute parameters
freqMin=94;
freqMax=102;
freqSteps=200;
%number of processors
nbProcs=20;
%parameters
%bounds of parameters
Ymin=1.5;
Ymax=2.75;
Rmin=0.1;
Rmax=1;
Xfixed=1.5;

for itN=[5 10 15 20 25 30]
    
    ns=itN;
    
    doeXmin=[Ymin Rmin];
    doeXmax=[Ymax Rmax];
    
    %
    mDOE=multiDOE(1,'IHS',ns,doeXmin,doeXmax);
    %
    samplePts=mDOE.sorted;
    %
    SILEX=wrapperSILEX;
    SILEX.resultFile='results_square/xfem_3_results.mat';
    SILEX.pythonCompute={'square_all.py'};
    SILEX.caseDefine='thick_x_up_wall';
    SILEX.freqMin=freqMin;
    SILEX.freqMax=freqMax;
    SILEX.nbSteps=freqSteps;
    SILEX.nbProc=nbProcs;
    %run on parameters
    SILEX.compute([Xfixed*ones(ns,1) samplePts]);
    %
    namefilsave=['results_square/' datestr(datetime,'YYYY-mm-DD_HH-MM-SS_') 'runDOE.mat'];
    save(namefilsave)
    
    itSILEX=1;
    %samplePts=SILEX.paraValFull;
    varResult=SILEX.varResult;
    resp=[];
    dresp=zeros(ns,2);
    FRF={};
    dFRF={};
    prefsquare=(20e-6)^2;
    for itS=1:size(samplePts,1)
        FRF{itS}=10*log10(varResult{itS}.AllFRF(2,:)/prefsquare);
        dFRFX=10*varResult{itS}.AllFRF(4,:)./varResult{itS}.AllFRF(2,:).*1/log(10);
        dFRFY=10*varResult{itS}.AllFRF(5,:)./varResult{itS}.AllFRF(2,:).*1/log(10);
        resp(itS)=mean(FRF{itS});
        dresp(itS,1)=mean(dFRFX);
        dresp(itS,2)=mean(dFRFY);
    end
    
    namefilsave=['results/' datestr(datetime,'YYYY-mm-DD_HH-MM-SS_') 'postProc.mat'];
    save(namefilsave)
    
    nP=30;
    XX=linspace(Ymin,Ymax,nP);
    YY=linspace(Rmin,Rmax,nP);
    [XXm,YYm]=meshgrid(XX,YY);
    gridRef=[];
    gridRef(:,:,1)=XXm;
    gridRef(:,:,2)=YYm;
    
    
    Meta=GRENAT('GKRG',samplePts,resp(:),dresp);
    Meta.train;
    
    MetaB=GRENAT('KRG',samplePts,resp(:));
    MetaB.train;
    
    Meta.defineRef(gridRef);
    Meta.eval(gridRef);
    Meta.confDisp.newFig=false;
    Meta.show;
    Meta.confDisp.conf('samplePts',true);
    figure;Meta.showResp;
    f=xlabel('$y_w$ [m]','Interpreter','latex')
    ylabel('$R_w$ [m]','Interpreter','latex')
    zlabel('$\tilde{L_p}$ [dB]','Interpreter','latex')
    title('CKRG')
    %view([63 52])
    hlight=light;               % active light
    lighting('gouraud');         % type of rendering
    lightangle(hlight,48,70);    % direction of the light
    saveas(f,['WCCM/CKRG_' num2str(ns) '_94-102.eps'],'epsc')
    
    
    MetaB.defineRef(gridRef);
    MetaB.eval(gridRef);
    MetaB.confDisp.newFig=false;
    MetaB.show;
    MetaB.confDisp.conf('samplePts',true);
    f=figure;MetaB.showResp;
    xlabel('$y_w$ [m]','Interpreter','latex')
    ylabel('$R_w$ [m]','Interpreter','latex')
    zlabel('$\tilde{L_p}$ [dB]','Interpreter','latex')
    %view([63 52])
    hlight=light;               % active light
    lighting('gouraud');         % type of rendering
    lightangle(hlight,48,70);    % direction of the light
    title('KRG')
    saveas(f,['WCCM/KRG_' num2str(ns) '_94-102.eps'],'epsc')
    
    %load ref
    REF=load('2018-07-19_15-21-17_xfem_3_Xfixed_94_102_900');
    %regroup data
    funC=@(X) X.AllFRF(2,:);
    tmpV=cellfun(funC,REF.varResult,'UniformOutput',false);
    FrFallRAW=vertcat(tmpV{:});
    %gradients
    funC=@(X) X.AllFRF(3,:);
    tmpV=cellfun(funC,REF.varResult,'UniformOutput',false);
    FrFallG1RAW=vertcat(tmpV{:});
    funC=@(X) X.AllFRF(4,:);
    tmpV=cellfun(funC,REF.varResult,'UniformOutput',false);
    FrFallG2RAW=vertcat(tmpV{:});
    %find NaN
    maskNAN=isnan(FrFallRAW);
    paraNAN=REF.paraValFull(maskNAN(:,end),:);
    if ~isempty(paraNAN);fprintf('NAN values for parameters: \n');end
    fprintf('%g %g\n',paraNAN')
    %remove NaN
    FrFallRAW(maskNAN(:,end),:)=[];
    % compute in dB
    prefsquare=(20e-6)^2;
    FrFall=10*log10(FrFallRAW./prefsquare);
    FrFallG1=10*FrFallG1RAW./FrFallRAW.*1/log(10);
    FrFallG2=10*FrFallG2RAW./FrFallRAW.*1/log(10);
    LpMean=zeros(REF.nP);
    LpMean(:)=mean(FrFall,2);
    
    
    gridBRef=[];
    gridBRef(:,:,1)=REF.Ym;
    gridBRef(:,:,2)=REF.Rm;
    
    
    Meta.defineRef('sampleRef',gridBRef,'respRef',LpMean);
    aerr=Meta.err;
    
    MetaB.defineRef('sampleRef',gridBRef,'respRef',LpMean);
    berr=MetaB.err;
    
    
    namefilsave=['results/' datestr(datetime,'YYYY-mm-DD_HH-MM-SS_') 'META.mat'];
    save(namefilsave)
    
end

%load ref
REF=load('2018-07-19_15-21-17_xfem_3_Xfixed_94_102_900');
%regroup data
funC=@(X) X.AllFRF(2,:);
tmpV=cellfun(funC,REF.varResult,'UniformOutput',false);
FrFallRAW=vertcat(tmpV{:});
%gradients
funC=@(X) X.AllFRF(3,:);
tmpV=cellfun(funC,REF.varResult,'UniformOutput',false);
FrFallG1RAW=vertcat(tmpV{:});
funC=@(X) X.AllFRF(4,:);
tmpV=cellfun(funC,REF.varResult,'UniformOutput',false);
FrFallG2RAW=vertcat(tmpV{:});
%find NaN
maskNAN=isnan(FrFallRAW);
paraNAN=REF.paraValFull(maskNAN(:,end),:);
if ~isempty(paraNAN);fprintf('NAN values for parameters: \n');end
fprintf('%g %g\n',paraNAN')
%remove NaN
FrFallRAW(maskNAN(:,end),:)=[];
% compute in dB
prefsquare=(20e-6)^2;
FrFall=10*log10(FrFallRAW./prefsquare);
FrFallG1=10*FrFallG1RAW./FrFallRAW.*1/log(10);
FrFallG2=10*FrFallG2RAW./FrFallRAW.*1/log(10);
LpMean=zeros(REF.nP);
LpMean(:)=mean(FrFall,2);
%%plot surfaces
f=figure
LpMean=zeros(REF.nP);
LpMean(:)=mean(FrFall,2);
surf(REF.Ym,REF.Rm,LpMean)
hlight=light;               % active light
lighting('gouraud');         % type of rendering
lightangle(hlight,48,70);    % direction of the light
xlabel('$y_w$ [m]','Interpreter','latex');
ylabel('$R_w$ [m]','Interpreter','latex');
zlabel('$\overline{L_p}$ [dB]','Interpreter','latex');
colorbar
%view([63 52])
saveas(f,'WCCM/Ref_obj_94-102.eps','epsc')
f=figure
LpMean=zeros(REF.nP);
LpMean(:)=mean(FrFallG1,2);
surf(REF.Ym,REF.Rm,LpMean)
hlight=light;               % active light
lighting('gouraud');         % type of rendering
lightangle(hlight,48,70);    % direction of the light
xlabel('$y_w$ [m]','Interpreter','latex');
ylabel('$R_w$ [m]','Interpreter','latex');
zlabel('$\frac{\partial\overline{L_p}}{\partial x_w}$ [dB/m]','Interpreter','latex');
colorbar
%view([63 52])
saveas(f,'WCCM/Ref_GX_94-102.eps','epsc')
f=figure;
LpMean=zeros(REF.nP);
LpMean(:)=mean(FrFallG2,2);
surf(REF.Ym,REF.Rm,LpMean)
hlight=light;               % active light
lighting('gouraud');         % type of rendering
lightangle(hlight,48,70);    % direction of the light
xlabel('$y_w$ [m]','Interpreter','latex');
ylabel('$R_w$ [m]','Interpreter','latex');
zlabel('$\frac{\partial\overline{L_p}}{\partial y_w}$ [dB/m]','Interpreter','latex');
colorbar
%view([63 52])
saveas(f,'WCCM/Ref_GY_94-102.eps','epsc')
