%parametric study
clear all
%number of processors
nbProcs=20;
%frequency range
freq.max=100;
freq.min=90;
freq.steps=100;
%number of step per parameter 
nP=30;

%bounds of parameters
angle.min=0;
angle.max=180;
position.min=0.1;
position.max=2.7;

%%%%
[XX,YY]=meshgrid(linspace(angle.min,angle.max,nP),linspace(position.min,position.max,nP));
paraVal(:,:,1)=XX;
paraVal(:,:,2)=YY;
%%%%
%start wrapper
SILEX=wrapperSILEX;
SILEX.nbSteps=freq.steps;
SILEX.freqMax=freq.max;
SILEX.freqMin=freq.min;
SILEX.nbProc=nbProcs;
%run on parameters
SILEX.compute(paraVal);
%save everything
save(SILEX.saveFileFull,'-append');
