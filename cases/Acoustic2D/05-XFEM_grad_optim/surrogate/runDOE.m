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

%ii=1;
for ii=1:numel(nscomb)
    ns=nscomb(ii);%40;
    for jj=1:nbRep        
        %
        mDOE=multiDOE(2,'LHS',ns,Xmin,Xmax);
        %
        samplePts{ii,jj}=mDOE.sorted;
        %
        SILEX(ii,jj)=wrapperSILEX;
        SILEX(ii,jj).freqMin=freqMin;
        SILEX(ii,jj).freqMax=freqMax;
        SILEX(ii,jj).nbSteps=freqSteps;
        SILEX(ii,jj).nbProc=nbProcs;
        %run on parameters
        SILEX(ii,jj).compute(samplePts{ii,jj});
    %
    namefilsave=[datestr(datetime,'YYYY-mm-DD_HH-MM-SS_') 'runDOE.mat'];
    save(namefilsave)
    end
end

namefilsave=[datestr(datetime,'YYYY-mm-DD_HH-MM-SS_') 'runDOE.mat'];
save(namefilsave)
