%% execute SILEX for building kriging surrogate model
clear all
%add path
addpath('grenat')

%load folder structure
initDirGRENAT;
%
nsmin=20;
nsmax=40;
nbRep=1;
nscomb=nsmin:nsmax;
%compute parameters
freqMin=90;
freqMax=100;
freqSteps=100;
%number of processors
nbProcs=20;
%parameters
angle.min=0;
angle.max=180;
position.min=0.1;
position.max=2.7;
Xmin=[angle.min position.min];
Xmax=[angle.max,position.max];

ns=20;

%
mDOE=multiDOE(2,'IHS',ns,Xmin,Xmax);
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
namefilsave=[datestr(datetime,'YYYY-mm-DD_HH-MM-SS_') 'runDOE.mat'];
save(namefilsave)
    
itSILEX=1;
samplePts=SILEX.paraValFull;
resp=[];
for itS=1:size(samplePts,1)
    prefsquare=(20e-6)^2;
    FRF=10*log10(SILEX.varResult{itS}.AllFRF(2,:)/prefsquare);
    resp(itS)=mean(FRF);
end

nP=30;
[XX,YY]=meshgrid(linspace(angle.min,angle.max,nP),linspace(position.min,position.max,nP));
gridRef(:,:,1)=XX;
gridRef(:,:,2)=YY;
   
addpath('grenat')

Meta=GRENAT('RBF',samplePts,resp(:));
Meta.train;

Meta.defineRef(gridRef);
Meta.eval(gridRef);
Meta.show;   
 

namefilsave=[datestr(datetime,'YYYY-mm-DD_HH-MM-SS_') 'META.mat'];
save(namefilsave)