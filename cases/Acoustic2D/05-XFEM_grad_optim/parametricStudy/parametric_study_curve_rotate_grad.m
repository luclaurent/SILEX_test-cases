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
    1.5 4 0.75
    1.5+decal 4 0.75
    ];

paraInit=[1.5 4 0.75
    1.5+decal 4  0.75
    1.5 4+decal 0.75
    1.5 4 0.75+decal
    1.5-decal 4  0.75
    1.5 4-decal 0.75
    1.5 4 0.75-decal
    ];

%start wrapper
SILEX=wrapperSILEXcurve_rotate;
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

if v==1
    FRFX=varResult{2}.AllFRF(2,:);
    FRFY=varResult{3}.AllFRF(2,:);
    FRFR=varResult{4}.AllFRF(2,:);
    %
    dFRFXFD=(FRFX-FRFI)./decal;
    dFRFYFD=(FRFY-FRFI)./decal;
    dFRFRFD=(FRFR-FRFI)./decal;
elseif v==2
    FRFI=varResult{1}.AllFRF(2,:);
    FRFXp=varResult{2}.AllFRF(2,:);
    FRFYp=varResult{3}.AllFRF(2,:);
    FRFRp=varResult{4}.AllFRF(2,:);
    FRFXm=varResult{5}.AllFRF(2,:);
    FRFYm=varResult{6}.AllFRF(2,:);
    FRFRm=varResult{7}.AllFRF(2,:);
    %
    dFRFXFD=(FRFXp-FRFXm)./(2*decal);
    dFRFYFD=(FRFYp-FRFYm)./(2*decal);
    dFRFRFD=(FRFRp-FRFRm)./(2*decal);
end


figure;
plot(freq,FRFI)
hold on
plot(freq,FRFX)
plot(freq,FRFY)
plot(freq,FRFR)
figure;
plot(freq,12.53*dFRFX);hold on
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
