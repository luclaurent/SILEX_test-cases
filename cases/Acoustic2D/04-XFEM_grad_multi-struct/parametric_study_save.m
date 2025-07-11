%parametric study
clear all
close all
%number of processors
nbProcs=20;
%frequency range
freq.max=1000;%300;
freq.min=10;%10;
freq.steps=1000;
%number of step per parameter 
nP=2;
%bounds of parameters
position.min=0.505;
position.max=0.695;

%%%%
XX=linspace(position.min,position.max,nP);
XXb=XX+1e-5;
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
%start wrapper
SILEXb=wrapperSILEX;
SILEXb.nbSteps=freq.steps;
SILEXb.freqMax=freq.max;
SILEXb.freqMin=freq.min;
SILEXb.nbProc=nbProcs;
%run on parameters
SILEXb.compute(XXb');
save(SILEXb.saveFileFull,'-append')

prefsquare=(20e-6)^2;
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
sPtsb=SILEXb.paraValFull;
sPts=SILEX.paraValFull;
lF=SILEX.varResult{1}.AllFRF(1,:);
figure
listPlot=1;%1:size(sPts,1);
for itS=listPlot
    FRF=10*log10(SILEX.varResult{itS}.AllFRF(2,:)./prefsquare);
    FRFb=10*log10(SILEXb.varResult{itS}.AllFRF(2,:)./prefsquare);
    dFRF=10*SILEX.varResult{itS}.AllFRF(3,:)./SILEX.varResult{itS}.AllFRF(2,:).*1/log(10);
    dFRFFD=(FRFb-FRF)/(sPtsb(itS)-sPts(itS));
    subplot(411)
    plot(lF,FRF)
    title('FRF')
    hold on
    subplot(412)
    plot(lF,dFRF)
    title('dFRF')
    hold on
    subplot(413)
    plot(lF,dFRFFD)
    title('dFRF FD')
    hold on
    subplot(414)
    plot(lF,(dFRF-dFRFFD)./dFRFFD)
    title('diff');
    hold on
end
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%




samplePts=SILEX.paraValFull;
prefsquare=(20e-6)^2;

%%search for larger magnitude on objective function
%lower/larger frequencies
lF=SILEX.varResult{1}.AllFRF(1,:);
loF=min(lF);
hiF=max(lF);
%step frequency
sF=lF(2)-loF;
%loop for window
winMin=10;
deltaFRF=zeros(numel(lF)-winMin);
loFTable=deltaFRF;
hiFTable=deltaFRF;
if false
for itL=1:numel(lF)-winMin
    for itH=(itL+winMin):numel(lF)
        loFtmp=lF(itL);
        hiFtmp=lF(itH);
        fprintf('lo %g hi %g\n',loFtmp,hiFtmp);
        for itS=1:size(samplePts,1)
            FRF=10*log10(SILEX.varResult{itS}.AllFRF(2,itL:itH)./prefsquare);
            valObj(itS)=mean(FRF);
        end
        deltaFRF(itL,itH)=max(valObj)-min(valObj);
        loFTable(itL,itH)=loFtmp;
        hiFTable(itL,itH)=hiFtmp;
    end
end
end
%max of deltaFRF
[objMax,Imax]=max(deltaFRF(:));
loFMax=loFTable(Imax);
hiFMax=hiFTable(Imax);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot FRF in window
loFW=645;%645;%86;%449;%239;
hiFW=767;%767;%255;%543;%386;
%find indexes
[loFexact,Ilo]=min(abs(lF-loFW));
[hiFexact,Ihi]=min(abs(lF-hiFW));
lfW=lF(Ilo:Ihi);
%compute FRF
samplePtsPlot=samplePts(:);%:5);
valObjW=zeros(1,size(samplePtsPlot,1));
valDObjW=valObjW;
figure
for itS=1:size(samplePtsPlot,1)
    FRF=10*log10(SILEX.varResult{itS}.AllFRF(2,Ilo:Ihi)./prefsquare);
    dFRF=10*SILEX.varResult{itS}.AllFRF(3,Ilo:Ihi)./SILEX.varResult{itS}.AllFRF(2,Ilo:Ihi).*1/log(10);
    %FRF=SILEX.varResult{itS}.AllFRF(2,Ilo:Ihi);
    %dFRF=SILEX.varResult{itS}.AllFRF(3,Ilo:Ihi);
    %dFRFFD=(SILEX.varResult{2}.AllFRF(2,Ilo:Ihi)-SILEX.varResult{1}.AllFRF(2,Ilo:Ihi))/1e-5;
    valObjW(itS)=mean(FRF);
    valDObjW(itS)=mean(dFRF);
    %subplot(121)
    plot(lfW,(FRF-mean(FRF))./(max(FRF)-mean(FRF)),'r');
    hold on
    %subplot(122)
    plot(lfW,dFRF./max(abs(dFRF)),'k');
   % plot(lfW,dFRFFD./max(abs(dFRFFD)));
    line([min(lfW) max(lfW)],[0 0])
    hold on
end
figure;
plot(samplePtsPlot,valObjW);
figure
vD=(valObjW(2:end)-valObjW(1:end-1))/(samplePtsPlot(2)-samplePtsPlot(1));
plot(samplePtsPlot,valDObjW)
hold on
plot(samplePtsPlot,vD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



figure
plot(SILEX.varResult{1}.AllFRF(1,:),SILEX.varResult{1}.AllFRF(2,:),'-b')
hold on
plot(SILEX.varResult{1}.AllFRF(1,:),SILEX.varResult{1}.AllFRF(3,:),'-r')
hold on


figure
for itS=1:size(samplePts,1)
plot(SILEX.varResult{itS}.AllFRF(1,:),SILEX.varResult{itS}.AllFRF(2,:))
hold on
end
hold off
figure
for itS=1:size(samplePts,1)
plot(SILEX.varResult{itS}.AllFRF(1,:),SILEX.varResult{itS}.AllFRF(3,:),'-r')
hold on
end
hold off



samplePts=SILEX.paraValFull;
resp=[];
dresp=[];
FRF={};
dFRF={};
for itS=1:size(samplePts,1)
    FRF{itS}=10*log10(SILEX.varResult{itS}.AllFRF(2,:)./prefsquare);
    dFRF{itS}=SILEX.varResult{itS}.AllFRF(3,:)./log(10)./SILEX.varResult{itS}.AllFRF(2,:);
    resp(itS)=mean(FRF{itS});
    dresp(itS)=mean(dFRF{itS});
end

figure
for itS=1:size(samplePts,1)
plot(SILEX.varResult{itS}.AllFRF(1,:),FRF{itS})
hold on
end
hold off
figure
for itS=1:size(samplePts,1)
plot(SILEX.varResult{itS}.AllFRF(1,:),dFRF{itS})
hold on
end
hold off


figure;
subplot(1,2,1)
plot(samplePts,resp)
subplot(1,2,2)
plot(samplePts,dresp,'-r')