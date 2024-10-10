% fun acoustic

function [Z,GZ]=funSILEX_Rfixed(X,R)
Wmax=82;
Wmin=49;

Rs=repmat(R,size(X,1),1);
    

SILEX=wrapperSILEX;

SILEX.resultFile='results_square/xfem_3_results.mat';
SILEX.pythonCompute={'square_all.py'};
SILEX.caseDefine='thick_x_up_wall';
SILEX.freqMin=Wmin;
SILEX.freqMax=Wmax;
SILEX.nbSteps=200;
SILEX.nbProc=20;

SILEX.compute([X Rs]);


prefsquare=(20e-6)^2;

nbZ=numel(SILEX.varResultFinal);
Z=zeros(nbZ,1);
for it=1:nbZ
    FRF=10*log10(SILEX.varResultFinal{it}.AllFRF(2,:)/prefsquare);
    dFRFX=10*SILEX.varResultFinal{it}.AllFRF(3,:)./SILEX.varResultFinal{it}.AllFRF(2,:).*1/log(10);
    dFRFY=10*SILEX.varResultFinal{it}.AllFRF(4,:)./SILEX.varResultFinal{it}.AllFRF(2,:).*1/log(10);
    Z(it)=mean(FRF);
    GZ(it,1)=mean(dFRFX);
    GZ(it,2)=mean(dFRFY);
end
end