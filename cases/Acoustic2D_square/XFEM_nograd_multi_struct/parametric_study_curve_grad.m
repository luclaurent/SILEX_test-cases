%parametric study
clear all
close all
%number of processors
nbProcs=20;
%frequency range
freq.max=12;%300;
freq.min=10;%10;
freq.steps=500;
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
    4 1.5 0.75
    4 1.5+decal 0.75
    ];

paraInit=[
    4 1.5 0.75
    4+decal 1.5 0.75
    4 1.5+decal 0.75
    4 1.5 0.75+decal
    ];

%start wrapper
SILEX=wrapperSILEXcurve;
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
FRFX=varResult{2}.AllFRF(2,:);
FRFY=varResult{3}.AllFRF(2,:);
FRFR=varResult{4}.AllFRF(2,:);
%
dFRFX=varResult{1}.AllFRF(3,:);
dFRFY=varResult{1}.AllFRF(4,:);
dFRFR=varResult{1}.AllFRF(5,:);
%
dFRFXFD=(FRFX-FRFI)./decal;
dFRFYFD=(FRFY-FRFI)./decal;
dFRFRFD=(FRFR-FRFI)./decal;


figure;
plot(freq,FRFI)
hold on
plot(freq,FRFX)
plot(freq,FRFY)
plot(freq,FRFR)
figure;
plot(freq,dFRFX);hold on
plot(freq,dFRFXFD,'LineWidth',2);
legend('X','FD X')
figure;
plot(freq,dFRFY);hold on
plot(freq,dFRFYFD,'LineWidth',2);
legend('Y','FD Y')
figure;
plot(freq,dFRFR);hold on
plot(freq,dFRFRFD,'LineWidth',2);
legend('R','FD R')

figure;
plot(freq,dFRFX-dFRFXFD,'LineWidth',2);
legend('X')
figure;
plot(freq,dFRFY-dFRFYFD,'LineWidth',2);
legend('Y')
figure;
plot(freq,dFRFR-dFRFRFD,'LineWidth',2);
legend('R')
