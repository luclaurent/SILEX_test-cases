%parametric study
clear all
close all
%number of processors
nbProcs=20;
%frequency range
freq.max=20;%300;
freq.min=10;%10;
freq.steps=200;
%number of step per parameter 
nP=1;
%bounds of parameters
position.min=0.505*7;
position.max=0.695*7;

%%%%
decal=1e-5;
XX=linspace(position.min,position.max,nP);
XXb=XX+decal;
paraVal=[XX XXb]';

paraVal=[
    1.5 4 0.75
    1.5+decal 4 0.75
    ];

paraRAW=[2.5 2.5 0.1 10];

paraInit=repmat(paraRAW,9,1)+[
    0 0 0 0
    decal 0 0 0
    0 decal 0 0
    0 0 decal 0
    0 0 0 decal
    -decal 0 0 0 
    0 -decal 0 0
    0  0 -decal 0
    0 0 0 -decal
    ];

%start wrapper
SILEX=wrapperSILEX;
SILEX.resultFile='results_square/xfem_3_results.mat';
SILEX.pythonCompute={'square_u_custom.py'};
SILEX.nbSteps=freq.steps;
SILEX.freqMax=freq.max;
SILEX.freqMin=freq.min;
SILEX.nbProc=nbProcs;
%run on parameters
SILEX.compute(paraInit);


varResult=SILEX.varResult;
paraValFull=SILEX.paraValFull;
save(SILEX.saveFileFull,'-append')


freq=varResult{1}.AllFRF(1,:);
FRFI=varResult{1}.AllFRF(2,:);
dFRFX=varResult{1}.AllFRF(3,:);
dFRFY=varResult{1}.AllFRF(4,:);
dFRFR=varResult{1}.AllFRF(5,:);
dFRFT=varResult{1}.AllFRF(6,:);

v=2;


if v==1
    FRFX=varResult{2}.AllFRF(2,:);
    FRFY=varResult{3}.AllFRF(2,:);
    FRFR=varResult{4}.AllFRF(2,:);
    %
    dFRFXFD=(FRFX-FRFI)./decal;
    dFRFYFD=(FRFY-FRFI)./decal;
    dFRFRFD=(FRFR-FRFI)./decal;
elseif v==2
    FRFXp=varResult{2}.AllFRF(2,:);
    FRFX=FRFXp;
    FRFYp=varResult{3}.AllFRF(2,:);
    FRFY=FRFYp;
    FRFRp=varResult{4}.AllFRF(2,:);
    FRFR=FRFRp;
    FRFTp=varResult{5}.AllFRF(2,:);
    FRFT=FRFTp;
    FRFXm=varResult{6}.AllFRF(2,:);
    FRFYm=varResult{7}.AllFRF(2,:);
    FRFRm=varResult{8}.AllFRF(2,:);
    FRFTm=varResult{9}.AllFRF(2,:);
    %
    dFRFXFD=(FRFXp-FRFXm)./(2*decal);
    dFRFYFD=(FRFYp-FRFYm)./(2*decal);
    dFRFRFD=(FRFRp-FRFRm)./(2*decal);
    dFRFTFD=(FRFTp-FRFTm)./(2*decal);
    %
    dFRFXFDa=(FRFXp-FRFI)./decal;
    dFRFYFDa=(FRFYp-FRFI)./decal;
    dFRFRFDa=(FRFRp-FRFI)./decal;
    dFRFTFDa=(FRFTp-FRFI)./decal;
end


figure;
plot(freq,FRFI)
hold on
plot(freq,FRFXp)
plot(freq,FRFYp)
plot(freq,FRFRp)
plot(freq,FRFTp)
figure;
plot(freq,dFRFX);hold on
plot(freq,dFRFXFD,'LineWidth',2);
plot(freq,dFRFXFDa,'LineWidth',2);
legend('X','FD X','FD1 X')
figure;
plot(freq,dFRFY);hold on
plot(freq,dFRFYFD,'LineWidth',2);
plot(freq,dFRFYFDa,'LineWidth',2);
legend('Y','FD Y','FD1 Y')
figure;
plot(freq,dFRFR);hold on
plot(freq,dFRFRFD,'LineWidth',2);
plot(freq,dFRFRFDa,'LineWidth',2);
legend('R','FD R','FD1 R')
figure;
plot(freq,dFRFT);hold on
plot(freq,dFRFTFD,'LineWidth',2);
plot(freq,dFRFTFDa,'LineWidth',2);
legend('T','FD T','FD1 T')
figure;
plot(freq,dFRFX);
hold on
plot(freq,dFRFY);
plot(freq,dFRFR);
plot(freq,dFRFT);
legend('X','Y','R','T')
figure;
plot(freq,dFRFXFD);
hold on
plot(freq,dFRFYFD);
plot(freq,dFRFRFD);
plot(freq,dFRFTFD);
legend('FD X','FD Y','FD R','FD T')

% figure;
% plot(freq,dFRFX-dFRFXFD,'LineWidth',2);
% legend('X')
% figure;
% plot(freq,dFRFY-dFRFYFD,'LineWidth',2);
% legend('Y')
% figure;
% plot(freq,dFRFR-dFRFRFD,'LineWidth',2);
% legend('R')