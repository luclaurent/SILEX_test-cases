%parametric study
clear all
close all
%number of processors
nbProcs=20;
%frequency range
freq.max=300;%300;
freq.min=10;%10;
freq.steps=500;
%number of step per parameter 
nP=2;
%bounds of parameters
position.min=0.505*7;
position.max=0.695*7;

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



samplePts=SILEX.paraValFull;
prefsquare=(20e-6)^2;

%%search for larger magnitude on objective function
%lower/larger frequencies
lF=SILEX.varResult{1}{1}.AllFRF(1,:);
loF=min(lF);
hiF=max(lF);
%step frequency
sF=lF(2)-loF;
%loop for window
winMin=10;
deltaFRF=zeros(numel(lF)-winMin);
loFTable=deltaFRF;
hiFTable=deltaFRF;
if true%false
for itL=1:numel(lF)-winMin
    for itH=(itL+winMin):numel(lF)
        loFtmp=lF(itL);
        hiFtmp=lF(itH);
        fprintf('lo %g hi %g\n',loFtmp,hiFtmp);
        for itS=1:size(samplePts,1)
            FRF=10*log10(SILEX.varResult{1}{itS}.AllFRF(2,itL:itH)./prefsquare);
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
target=10;
nbTest=10;
[~,Is]=sort(abs(deltaFRF(:)-target));
valObj=[];
figure
for ii=1:nbTest
    loFtmp=loFTable(Is(ii));
    hiFtmp=hiFTable(Is(ii));
    %find indexes
    [~,Ilo]=min(abs(lF-loFtmp));
    [~,Ihi]=min(abs(lF-hiFtmp));
    loFexact=lF(Ilo);
    hiFexact=lF(Ihi);
    for itS=1:size(samplePts,1)
        FRF=10*log10(SILEX.varResult{1}{itS}.AllFRF(2,Ilo:Ihi)./prefsquare);
        valObj(itS)=mean(FRF);
    end
    plot(samplePts,valObj,'DisplayName',['lo ' num2str(loFexact) 'hi ' num2str(hiFexact)])
    hold on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot FRF in window
allVal=[
    77 90
    122 134
    177 187
    194 204
    179 190
%     662 888
%     507 678
%     329 679
%     401 585
%     413 596
%     589 816
%     529 744
%     662 888
%     677 906
%     668 812
%     705 857
%     788 832
%     392 444
%     732 774
%     692 729
%     444 469
    ];
for ii=1:5%16
    valC=ii;
    loFW=allVal(valC,1);%507;%645;%86;%449;%239;
    hiFW=allVal(valC,2);%888;%767;%255;%543;%386;
    %find indexes
    [loFexact,Ilo]=min(abs(lF-loFW));
    [hiFexact,Ihi]=min(abs(lF-hiFW));
    loFexact=lF(Ilo);
    hiFexact=lF(Ihi);
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
    title(['lo ' num2str(loFexact) 'hi ' num2str(hiFexact)]);
    figure
    vD=(valObjW(2:end)-valObjW(1:end-1))/(samplePtsPlot(2)-samplePtsPlot(1));
    plot(samplePtsPlot,valDObjW)
    hold on
    plot(samplePtsPlot(2:end),vD)
end
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