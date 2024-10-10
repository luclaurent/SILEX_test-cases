%% execute SILEX for building kriging surrogate model
clear all
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
freqMin=49;
freqMax=82;
freqSteps=200;
%number of processors
nbProcs=20;
%parameters
%bounds of parameters
Xmin=1;
Xmax=3.5;
Ymin=1.5;
Ymax=3.5;
Rmin=0.1;
Rmax=1;
Rfixed=0.5;

ns=10;

doeXmin=[Xmin Ymin];
doeXmax=[Xmax Ymax];

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
SILEX.compute([samplePts Rfixed*ones(ns,1)]);
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
    dFRFX=10*varResult{itS}.AllFRF(3,:)./varResult{itS}.AllFRF(2,:).*1/log(10);
    dFRFY=10*varResult{itS}.AllFRF(4,:)./varResult{itS}.AllFRF(2,:).*1/log(10);
    resp(itS)=mean(FRF{itS});
    dresp(itS,1)=mean(dFRFX);
    dresp(itS,2)=mean(dFRFY);
end

namefilsave=['results/' datestr(datetime,'YYYY-mm-DD_HH-MM-SS_') 'postProc.mat'];
save(namefilsave)

nP=30;
XX=linspace(Xmin,Xmax,nP);
YY=linspace(Ymin,Ymax,nP);
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
f=figure;Meta.showResp; 
xlabel('$x_w$ [m]','Interpreter','latex')
ylabel('$y_w$ [m]','Interpreter','latex')
zlabel('$\tilde{L_p}$ [dB]','Interpreter','latex')
title('CKRG')
view([63 52])
hlight=light;               % active light
lighting('gouraud');         % type of rendering
lightangle(hlight,48,70);    % direction of the light
saveas(f,['WCCM/CKRG_' num2str(ns) '_49-82.eps'],'epsc')


MetaB.defineRef(gridRef);
MetaB.eval(gridRef);
MetaB.confDisp.newFig=false;
MetaB.show;
MetaB.confDisp.conf('samplePts',true);
f=figure;MetaB.showResp;  
xlabel('$x_w$ [m]','Interpreter','latex')
ylabel('$y_w$ [m]','Interpreter','latex')
zlabel('$\tilde{L_p}$ [dB]','Interpreter','latex')
view([63 52])
hlight=light;               % active light
lighting('gouraud');         % type of rendering
lightangle(hlight,48,70);    % direction of the light
title('KRG')
saveas(f,['WCCM/KRG_' num2str(ns) '_49-82.eps'],'epsc')

%load ref
REF=load('2018-07-17_19-07-06_xfem_3_Rfixed_49_82_900');
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
gridBRef(:,:,1)=REF.Xm;
gridBRef(:,:,2)=REF.Ym;


Meta.defineRef('sampleRef',gridBRef,'respRef',LpMean);
aerr=Meta.err;

MetaB.defineRef('sampleRef',gridBRef,'respRef',LpMean);
berr=MetaB.err;
 

namefilsave=['WCCM/' datestr(datetime,'YYYY-mm-DD_HH-MM-SS_') 'META.mat'];
save(namefilsave)