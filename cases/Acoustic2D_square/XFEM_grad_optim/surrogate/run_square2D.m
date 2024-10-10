%parametric study
clear all
close all
%number of processors
nbProcs=20;
%frequency range
freq.max=82;%300;
freq.min=49;%10;
freq.steps=200;
%number of step per parameter 
nP=30;

%bounds of parameters
Xmin=1;
Xmax=2.75;
Ymin=1.5;
Ymax=3.5;
Rmin=0.1;
Rmax=1;

%%%%
YY=linspace(Ymin,Ymax,nP);
RR=0.5*ones(1,nP*nP);%linspace(Rmin,Rmax,nP);
XX=linspace(Xmin,Ymax,nP);%1.5*ones(1,nP*nP);
%
[Xm,Ym]=meshgrid(XX,YY);
paraVal=[Xm(:) Ym(:) RR(:)];
%start wrapper
SILEX=wrapperSILEX;
SILEX.resultFile='results_square/xfem_3_results.mat';
SILEX.pythonCompute={'square_all.py'};
SILEX.caseDefine='thick_x_up_wall';
SILEX.nbSteps=freq.steps;
SILEX.freqMax=freq.max;
SILEX.freqMin=freq.min;
SILEX.nbProc=nbProcs;
%run on parameterssave(SILEX.saveFileFull,'-append')
SILEX.compute(paraVal);

varResult=SILEX.varResult;
paraValFull=SILEX.paraValFull;
save(SILEX.saveFileFull,'-append')


% figure;
% for itV=1:numel(varResult)
%     plot(varResult{itV}.AllFRF(1,:),varResult{itV}.AllFRF(3,:),'DisplayName',num2str(itV));
%     hold on
% end

prefsquare=(20e-6)^2;

%compute objective function and gradients
valObj=zeros(nP);
valObjY=zeros(nP);
valObjR=zeros(nP);

for itS=1:size(paraValFull,1)
    FRF=10*log10(varResult{itS}.AllFRF(2,:)./prefsquare);
    dFRFX=10*varResult{itS}.AllFRF(3,:)./varResult{itS}.AllFRF(2,:).*1/log(10);
    dFRFY=10*varResult{itS}.AllFRF(4,:)./varResult{itS}.AllFRF(2,:).*1/log(10);
    valObj(itS)=mean(FRF);
    valObjY(itS)=mean(dFRFX);
    valObjR(itS)=mean(dFRFY);
end

figure;
subplot(131)
surf(Xm,Ym,valObj);
subplot(132)
surf(Xm,Ym,valObjY);
subplot(133)
surf(Xm,Ym,valObjR);
