%parametric study
clear all
close all
%number of processors
nbProcs=20;
%frequency range
freq.max=177;%300;
freq.min=134;%10;
freq.steps=200;
%number of step per parameter 
nP=30;
%bounds of parameters
Ymin=1;
Ymax=2.75;
Rmin=0.1;
Rmax=1;

%%%%
YY=linspace(Ymin,Ymax,nP);
RR=linspace(Rmin,Rmax,nP);
X=1.5*ones(1,nP*nP);
%
[Ym,Rm]=meshgrid(YY,RR);
paraVal=[X(:) Ym(:),Rm(:)];
%start wrapper
SILEX=wrapperSILEX;
SILEX.resultFile='results_square/xfem_3_results.mat';
SILEX.pythonCompute={'square_all.py'};
SILEX.caseDefine='thick_x_up_wall';
SILEX.nbSteps=freq.steps;
SILEX.freqMax=freq.max;
SILEX.freqMin=freq.min;
SILEX.nbProc=nbProcs;
%run on parameters
SILEX.compute(paraVal);

varResult=SILEX.varResult;
paraValFull=SILEX.paraValFull;
save(SILEX.saveFileFull,'-append')

samplePts=paraValFull;
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
surf(Ym,Rm,valObj);
subplot(132)
surf(Ym,Rm,valObjY);
subplot(133)
surf(Ym,Rm,valObjR);

% figure;
% for itV=1:numel(varResult)
%     plot(varResult{itV}.AllFRF(1,:),varResult{itV}.AllFRF(3,:),'DisplayName',num2str(itV));
%     hold on
% end
% 
% 
% samplePts=paraValFull;
% prefsquare=(20e-6)^2;
% 
% %%search for larger magnitude on objective function
% %lower/larger frequencies
% lF=varResult{1}.AllFRF(1,:);
% loF=min(lF);
% hiF=max(lF);
% %step frequency
% sF=lF(2)-loF;
% %loop for window
% winMin=10;
% deltaFRF=zeros(numel(lF)-winMin);
% loFTable=deltaFRF;
% hiFTable=deltaFRF;
% valObj=[]
% if true%false
% for itL=1:numel(lF)-winMin
%     for itH=(itL+winMin):numel(lF)
%         loFtmp=lF(itL);
%         hiFtmp=lF(itH);
%         fprintf('lo %g hi %g\n',loFtmp,hiFtmp);
%         for itS=1:size(samplePts,1)
%             FRF=10*log10(varResult{itS}.AllFRF(2,itL:itH)./prefsquare);
%             valObj(itS)=mean(FRF);
%         end
%         deltaFRF(itL,itH)=max(valObj)-min(valObj);
%         loFTable(itL,itH)=loFtmp;
%         hiFTable(itL,itH)=hiFtmp;
%     end
% end
% end
% %max of deltaFRF
% [objMax,Imax]=max(deltaFRF(:));
% loFMax=loFTable(Imax);
% hiFMax=hiFTable(Imax);
% 
% save(SILEX.saveFileFull,'-append')
% %%
% target=20;
% nbTest=10;
% [~,Is]=sort(abs(deltaFRF(:)-target));
% valObj=zeros(size(Ym));
% valObjX=zeros(size(Ym));
% valObjY=zeros(size(Ym));
% figure
% for ii=1:nbTest
%     loFtmp=loFTable(Is(ii));
%     hiFtmp=hiFTable(Is(ii));
%     %find indexes
%     [~,Ilo]=min(abs(lF-loFtmp));
%     [~,Ihi]=min(abs(lF-hiFtmp));
%     loFexact=lF(Ilo);
%     hiFexact=lF(Ihi);
%     for itS=1:size(samplePts,1)
%         FRF=10*log10(varResult{itS}.AllFRF(2,Ilo:Ihi)./prefsquare);
%         dFRFX=10*varResult{itS}.AllFRF(3,Ilo:Ihi)./varResult{itS}.AllFRF(2,Ilo:Ihi).*1/log(10);
%         dFRFY=10*varResult{itS}.AllFRF(4,Ilo:Ihi)./varResult{itS}.AllFRF(2,Ilo:Ihi).*1/log(10);
%         valObj(itS)=mean(FRF);
%         valObjX(itS)=mean(dFRFX);
%         valObjY(itS)=mean(dFRFY);
%     end
%     subplot(131)
%     surf(Ym,Rm,valObj,'DisplayName',['lo ' num2str(loFexact) 'hi ' num2str(hiFexact)])
%     hold on
%     subplot(132)
%     surf(Ym,Rm,valObjX,'DisplayName',['lo ' num2str(loFexact) 'hi ' num2str(hiFexact)])
%     hold on
%     subplot(133)
%     surf(Ym,Rm,valObjY,'DisplayName',['lo ' num2str(loFexact) 'hi ' num2str(hiFexact)])
%     hold on
% end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %plot FRF in window
% % allVal=[
% %     77 90
% %     122 134
% %     177 187
% %     194 204
% %     179 190
% % %     662 888
% % %     507 678
% % %     329 679
% % %     401 585
% % %     413 596
% % %     589 816
% % %     529 744
% % %     662 888
% % %     677 906
% % %     668 812
% % %     705 857
% % %     788 832
% % %     392 444
% % %     732 774
% % %     692 729
% % %     444 469
% %     ];
% % 
% % for ii=1:5%16
% %     valC=ii;
% %     loFW=allVal(valC,1);%507;%645;%86;%449;%239;
% %     hiFW=allVal(valC,2);%888;%767;%255;%543;%386;
% %     %find indexes
% %     [loFexact,Ilo]=min(abs(lF-loFW));
% %     [hiFexact,Ihi]=min(abs(lF-hiFW));
% %     loFexact=lF(Ilo);
% %     hiFexact=lF(Ihi);
% %     lfW=lF(Ilo:Ihi);
% %     %compute FRF    
% %     samplePtsPlot=samplePts;%:5);
% %     valObjW=zeros(1,size(samplePtsPlot,1));
% %     valDObjW=valObjW;
% %     figure
% %     for itS=1:size(samplePtsPlot,1)
% %         FRF=10*log10(varResult{itS}.AllFRF(2,Ilo:Ihi)./prefsquare);
% %         dFRF=10*varResult{itS}.AllFRF(3,Ilo:Ihi)./varResult{itS}.AllFRF(2,Ilo:Ihi).*1/log(10);
% %         %FRF=SILEX.varResult{itS}.AllFRF(2,Ilo:Ihi);
% %         %dFRF=SILEX.varResult{itS}.AllFRF(3,Ilo:Ihi);
% %         %dFRFFD=(SILEX.varResult{2}.AllFRF(2,Ilo:Ihi)-SILEX.varResult{1}.AllFRF(2,Ilo:Ihi))/1e-5;
% %         valObjW(itS)=mean(FRF);
% %         valDObjW(itS)=mean(dFRF);
% %         %subplot(121)
% %         plot(lfW,(FRF-mean(FRF))./(max(FRF)-mean(FRF)),'r');
% %         hold on
% %         %subplot(122)
% %         plot(lfW,dFRF./max(abs(dFRF)),'k');
% %         % plot(lfW,dFRFFD./max(abs(dFRFFD)));
% %         line([min(lfW) max(lfW)],[0 0])
% %         hold on
% %     end
% %     figure;
% %     plot(samplePtsPlot,valObjW);
% %     title(['lo ' num2str(loFexact) 'hi ' num2str(hiFexact)]);
% %     figure
% %     vD=(valObjW(2:end)-valObjW(1:end-1))/(samplePtsPlot(2)-samplePtsPlot(1));
% %     plot(samplePtsPlot,valDObjW)
% %     hold on
% %     plot(samplePtsPlot(2:end),vD)
% % end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % 
% % 
% % figure
% % plot(varResult{1}.AllFRF(1,:),varResult{1}.AllFRF(2,:),'-b')
% % hold on
% % plot(varResult{1}.AllFRF(1,:),varResult{1}.AllFRF(3,:),'-r')
% % hold on
% % 
