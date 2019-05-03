% fun acoustic

function [Z,GZ]=funSILEX_ZRfixed(X,Y)
Wmax=103;
Wmin=113;

Z=0.;
R=1.;

Zs=repmat(Z,size(X));
Rs=repmat(R,size(X));
    

SILEX=wrapperSILEX;

SILEX.resultFile='results/cavity_acou3D_struc_3D_v3results.mat';
SILEX.pythonCompute={'Main_acou3D_struc3D_v3_grad.py'};
SILEX.gradCompute=[0,1];
SILEX.freqMin=Wmin;
SILEX.freqMax=Wmax;
SILEX.nbSteps=50;
SILEX.nbProc=16;

SILEX.compute([X Y Zs Rs]);


prefsquare=(20e-6)^2;

nbZ=numel(SILEX.varResultFinal);
Z=zeros(nbZ,1);
GZ=zeros(nbZ,2);
for it=1:nbZ
    FRF=10*log10(SILEX.varResultFinal{it}.AllFRF(2,:)/prefsquare);
    dFRFX=10*SILEX.varResultFinal{it}.AllFRF(3,:)./SILEX.varResultFinal{it}.AllFRF(2,:).*1/log(10);
    dFRFY=10*SILEX.varResultFinal{it}.AllFRF(4,:)./SILEX.varResultFinal{it}.AllFRF(2,:).*1/log(10);
    Z(it)=mean(FRF);
    GZ(it,1)=mean(dFRFX);
    GZ(it,2)=mean(dFRFY);
end
end
