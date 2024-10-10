clear all
close all

addpath('grenat');initDirGRENAT

nsmin=2;
nsmax=10;
nbRep=10;
nscomb=nsmin:nsmax;
%parameters
%bounds of parameters
Xmin= 0.;
Xmax= 4.;
Ymin=0.;
Ymax=3.;

nP=30;
XX=linspace(Xmin,Xmax,nP);
YY=linspace(Ymin,Ymax,nP);
[XXm,YYm]=meshgrid(XX,YY);
gridRef=[];
gridRef(:,:,1)=XXm;
gridRef(:,:,2)=YYm;

%load ref
REF=load('paraXY_2500');
%load reference
refLp=REF.vF;
gridBRef(:,:,1)=REF.gX;
gridBRef(:,:,2)=REF.gY;

MetaREF=GRENAT('KRG',[REF.Xm(:) REF.Ym(:)],refLp(:));
MetaREF.confMeta.estimOn=false;
MetaREF.confMeta.cvOn=false;
MetaREF.confMeta.lVal=[1e-1 1e-1];
MetaREF.train;
MetaREF.defineRef(gridRef);
MetaREF.eval(gridRef);
MetaREF.confDisp.newFig=false;
MetaREF.show;
MetaREF.confDisp.conf('samplePts',true);
f=figure;MetaREF.showResp;

       %

Zref=zeros(size(XXm));
Zref(:)=MetaREF.eval([XXm(:) YYm(:)]);

for itR=1:nbRep
    for itN=[5:5:100]
        
        ns=itN;
        
        doeXmin=[Xmin Ymin];
        doeXmax=[Xmax Ymax];
        
        %
        mDOE=multiDOE(1,'IHS',ns,doeXmin,doeXmax);
        samplePts=mDOE.sorted;
        
        fileData=['data_' num2str(itN,'%03i') '_' num2str(itR,'%03i') '.mat'];
           
               
        folderSave='WCSMO';
        varX=('X');
        varY=('Y');
        freqMin=127;
        freqMax=137;
        
       
        [resp,dresp]=MetaREF.eval(samplePts);        
        
        Meta=GRENAT('GKRG',samplePts,resp(:),dresp);
        Meta.train;
        
        MetaB=GRENAT('KRG',samplePts,resp(:));
        MetaB.train;
        
        Meta.defineRef('sampleRef',gridRef,'respRef',Zref);
        Meta.eval(gridRef);
        Meta.confDisp.newFig=false;
        Meta.show;
        Meta.confDisp.conf('samplePts',true);
        f=figure;Meta.showResp;
        xlabel(['$' varX '$ [m]'],'Interpreter','latex')
        ylabel(['$' varY '$ [m]'],'Interpreter','latex')
        zlabel('$\tilde{L_p}$ [dB]','Interpreter','latex')
        title('CKRG')
        %view([63 52])
        hlight=light;               % active light
        lighting('gouraud');         % type of rendering
        lightangle(hlight,48,70);    % direction of the light
        saveas(f,fullfile(folderSave,['CKRG_' num2str(itN) '_' num2str(itR,'%03i') '_' num2str(freqMin) '-' num2str(freqMax) '.eps']),'epsc')
        
        
        MetaB.defineRef(gridRef);
        MetaB.eval(gridRef);
        MetaB.confDisp.newFig=false;
        MetaB.show;
        MetaB.confDisp.conf('samplePts',true);
        f=figure;MetaB.showResp;
        xlabel(['$' varX '$ [m]'],'Interpreter','latex')
        ylabel(['$' varY '$ [m]'],'Interpreter','latex')
        zlabel('$\tilde{L_p}$ [dB]','Interpreter','latex')
        %view([63 52])
        hlight=light;               % active light
        lighting('gouraud');         % type of rendering
        lightangle(hlight,48,70);    % direction of the light
        title('KRG')
        saveas(f,fullfile(folderSave,['KRG_' num2str(itN) '_' num2str(itR,'%03i') '_' num2str(freqMin) '-' num2str(freqMax) '.eps']),'epsc')
        
        
        
        % %regroup data
        % funC=@(X) X.AllFRF(2,:);
        % tmpV=cellfun(funC,REF.varResult,'UniformOutput',false);
        % FrFallRAW=vertcat(tmpV{:});
        % %gradients
        % funC=@(X) X.AllFRF(3,:);
        % tmpV=cellfun(funC,REF.varResult,'UniformOutput',false);
        % FrFallG1RAW=vertcat(tmpV{:});
        % funC=@(X) X.AllFRF(4,:);
        % tmpV=cellfun(funC,REF.varResult,'UniformOutput',false);
        % FrFallG2RAW=vertcat(tmpV{:});
        % %find NaN
        % maskNAN=isnan(FrFallRAW);
        % paraNAN=REF.paraValFull(maskNAN(:,end),:);
        % if ~isempty(paraNAN);fprintf('NAN values for parameters: \n');end
        % fprintf('%g %g\n',paraNAN')
        % %remove NaN
        % FrFallRAW(maskNAN(:,end),:)=[];
        % % compute in dB
        % prefsquare=(20e-6)^2;
        % FrFall=10*log10(FrFallRAW./prefsquare);
        % FrFallG1=10*FrFallG1RAW./FrFallRAW.*1/log(10);
        % FrFallG2=10*FrFallG2RAW./FrFallRAW.*1/log(10);
        % LpMean=zeros(REF.nP);
        % LpMean(:)=mean(FrFall,2);
        
        
        
        
        aerr=Meta.err;
        
        berr=MetaB.err;
        
        
        namefilsave=fullfile('results_surrogate',['META' fileData(5:end)]);
        save(namefilsave)
        
        close all
        
    end
end

%%plot surfaces
f=figure;
surf(REF.Xm,REF.Ym,refLp)
hlight=light;               % active light
lighting('gouraud');         % type of rendering
lightangle(hlight,48,70);    % direction of the light
xlabel(['$' varX '$ [m]'],'Interpreter','latex');
ylabel(['$' varY '$ [m]'],'Interpreter','latex');
zlabel('$\overline{L_p}$ [dB]','Interpreter','latex');
colorbar
%view([63 52])
saveas(f,fullfile(folderSave,['Ref_obj_' num2str(freqMin) '-' num2str(freqMax) '.eps']),'epsc')
f=figure;
surf(REF.Xm,REF.Ym,REF.gX)
hlight=light;               % active light
lighting('gouraud');         % type of rendering
lightangle(hlight,48,70);    % direction of the light
xlabel(['$' varX '$ [m]'],'Interpreter','latex');
ylabel(['$' varY '$ [m]'],'Interpreter','latex');
zlabel(['$\frac{\partial\overline{L_p}}{\partial ' varX '}$ [dB/m]'],'Interpreter','latex');
colorbar
%view([63 52])
saveas(f,fullfile(folderSave,['Ref_GX_' num2str(freqMin) '-' num2str(freqMax) '.eps']),'epsc')
f=figure;
surf(REF.Xm,REF.Ym,REF.gY)
hlight=light;               % active light
lighting('gouraud');         % type of rendering
lightangle(hlight,48,70);    % direction of the light
xlabel(['$' varX '$ [m]'],'Interpreter','latex');
ylabel(['$' varY '$ [m]'],'Interpreter','latex');
zlabel(['$\frac{\partial\overline{L_p}}{\partial ' varY '}$ [dB/m]'],'Interpreter','latex');
colorbar
%view([63 52])
saveas(f,fullfile(folderSave,['Ref_GY_' num2str(freqMin) '-' num2str(freqMax) '.eps']),'epsc')