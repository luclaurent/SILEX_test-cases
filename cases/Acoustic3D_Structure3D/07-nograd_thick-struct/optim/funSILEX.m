% fun acoustic

function Z=funSILEX(X)
Wmax=60;
Wmin=40;

classSILEX=wrapperSILEX;
classSILEX.compute(X);
nbZ=numel(classSILEX.varResultFinal);
Z=zeros(nbZ,1);
for it=1:nbZ
    listFreq=classSILEX.varResultFinal{it}.AllFRF(1,:);
    iXFreqW=find(listFreq>=Wmin&listFreq<=Wmax);
    Z(it)=mean(classSILEX.varResultFinal{it}.AllFRF(2,iXFreqW));
end
end