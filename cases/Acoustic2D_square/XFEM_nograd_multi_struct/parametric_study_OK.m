%parametric study
clear all
close all
%number of processors
nbProcs=20;
%frequency range
freq.max=10;%300;
freq.min=300;%10;
freq.steps=1000;
%number of step per parameter 
nP=50;
%bounds of parameters
position.min=0.505*7;
position.max=0.695*7;

%%%%
XX=linspace(position.min,position.max,nP);
paraVal=XX';
%start wrapper
SILEX=wrapperSILEX;
SILEX.nbSteps=freq.steps;
SILEX.freqMax=freq.max;
SILEX.freqMin=freq.min;
SILEX.nbProc=nbProcs;
%run on parameters
SILEX.compute(paraVal);
save(SILEX.saveFileFull,'-append')



samplePts=SILEX.paraValFull;
varResult=SILEX.varResult;
prefsquare=(20e-6)^2;

resp=[];
dresp=[];
FRF={};
dFRF={};
for itS=1:size(samplePts,1)
    FRF{itS}=10*log10(varResult{itS}.AllFRF(2,:)./prefsquare);
    dFRF{itS}=varResult{itS}.AllFRF(3,:)./log(10)./varResult{itS}.AllFRF(2,:);
    resp(itS)=mean(FRF{itS});
    dresp(itS)=mean(dFRF{itS});
end

figure
for itS=1:size(samplePts,1)
plot(varResult{itS}.AllFRF(1,:),FRF{itS})
hold on
end
hold off
figure
for itS=1:size(samplePts,1)
plot(varResult{itS}.AllFRF(1,:),dFRF{itS})
hold on
end
hold off


figure;
subplot(1,2,1)
plot(samplePts,resp)
subplot(1,2,2)
plot(samplePts,dresp,'-r')

save(SILEX.saveFileFull,'-append')