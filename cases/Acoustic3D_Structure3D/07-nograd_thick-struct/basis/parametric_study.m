%parametric study

%number of step per parameter 
nP=2;
%bounds of parameters
angle.min=0;
angle.max=180;
position.min=0.1;
position.max=2.7;

%%%%
[XX,YY]=meshgrid(linspace(angle.min,angle.max,nP),linspace(position.min,position.max,nP));
paraVal(:,:,1)=XX;
paraVal(:,:,2)=YY;
%%%%
%start wrapper
SILEX=wrapperSILEX;
%run on parameters
SILEX.compute(paraVal);

