%% execute SILEX for building kriging surrogate model
clear all
%add path
addpath('grenat')

%load folder structure
initDirGRENAT;
%
nsmin=10;
nsmax=40;
nscomb=nsmin:nsmax;
%parameters
angle.min=0;
angle.max=180;
position.min=0.1;
position.max=2.7;
Xmin=[angle.min position.min];
Xmax=[angle.max,position.max];
ii=1;
%for ii=1:numel(nscomb)
    ns=40;
    %
    mDOE=multiDOE(2,'IHS',ns,Xmin,Xmax);
    %
    samplePts{ii}=mDOE.sorted;
    %
    SILEX(ii)=wrapperSILEX;
    %run on parameters
    SILEX(ii).compute(samplePts{ii});
    %
%end

namefilsave=[datestr(datetime,'YYYY-mm-DD_HH-MM-SS_') 'runDOE.mat'];
save(namefilsave)