% Example of use of GMetaOpti
% L. LAURENT -- 04/09/2016 -- luc.laurent@lecnam.net

clear all
%add path
addpath('gmetaopti')

%load folder structure
initDirGMetaOpti;
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


%
%load the test function
testFun=@(X)funSILEX(X);
%testFun='funPeaks';
%type of optimization
typeOpti='EGO';
%init sample points
initSampleNs=15;

%load the optimisation object
optiO=GMetaOpti(typeOpti,Xmin,Xmax,testFun);
optiO.samplingNs=initSampleNs;
optiO.optiEGO;

namefilsave=[datestr(datetime,'YYYY-mm-DD_HH-MM-SS_') 'optiEGO.mat'];
save(namefilsave)