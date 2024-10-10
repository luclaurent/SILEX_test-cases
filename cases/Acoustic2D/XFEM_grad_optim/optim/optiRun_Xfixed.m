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
% Xmin=1;
% Xmax=3.5;
% Ymin=1.5;
% Ymax=3.5;
% Rmin=0.1;
% Rmax=1;
% Rfixed=0.5;

Ymin=1.5;
Ymax=2.75;
Rmin=0.1;
Rmax=1;
Xfixed=1.5;

doeXmin=[Ymin Rmin];
doeXmax=[Ymax Rmax];


%
%load the test function
testFun=@(X)funSILEX_Xfixed(X,Xfixed);
%testFun='funPeaks';
%type of optimization
typeOpti='EGO';
%init sample points
nSI=[5 10 15 20];% 30 35];
for iTs=1:numel(nSI)
    initSampleNs=nSI(iTs);
    
    %load the optimisation object
    optiO{iTs}=GMetaOpti(typeOpti,doeXmin,doeXmax,testFun);
    optiO{iTs}.metaType='GKRG';
    optiO{iTs}.optiGrad=true;
    optiO{iTs}.samplingNs=initSampleNs;
    optiO{iTs}.optiEGO;
end
namefilsave=[datestr(datetime,'YYYY-mm-DD_HH-MM-SS_') 'optiEGO.mat'];
save(namefilsave)
