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
freqMin=662;
freqMax=888;
freqSteps=200;
%number of processors
nbProcs=20;
%parameters
position.min=0.505;
position.max=0.695;
Xmin=position.min;
Xmax=position.max;

ns=5;

%
mDOE=multiDOE(1,'IHS',ns,Xmin,Xmax);
%
samplePts=mDOE.sorted;
%
SILEX=wrapperSILEX;
SILEX.freqMin=freqMin;
SILEX.freqMax=freqMax;
SILEX.nbSteps=freqSteps;
SILEX.nbProc=nbProcs;
%run on parameters
SILEX.compute(samplePts);
%
namefilsave=['results/' datestr(datetime,'YYYY-mm-DD_HH-MM-SS_') 'runDOE.mat'];
save(namefilsave)
    
itSILEX=1;
samplePts=SILEX.paraValFull;
varResult=SILEX.varResult;
resp=[];
dresp=[];
FRF={};
dFRF={};
prefsquare=(20e-6)^2;
for itS=1:size(samplePts,1)
    FRF{itS}=10*log10(varResult{itS}.AllFRF(2,:)/prefsquare);
    dFRF{itS}=varResult{itS}.AllFRF(3,:)./log(10)./varResult{itS}.AllFRF(2,:);
    resp(itS)=mean(FRF{itS});
    dresp(itS)=mean(dFRF{itS});
end

namefilsave=['results/' datestr(datetime,'YYYY-mm-DD_HH-MM-SS_') 'postProc.mat'];
save(namefilsave)

nP=100;
XX=linspace(position.min,position.max,nP);
gridRef(:,:,1)=XX;
   

Meta=GRENAT('GKRG',samplePts,resp(:),dresp(:));
Meta.train;

Meta.defineRef(gridRef);
Meta.eval(gridRef);
Meta.show;   
 

namefilsave=['results/' datestr(datetime,'YYYY-mm-DD_HH-MM-SS_') 'META.mat'];
save(namefilsave)