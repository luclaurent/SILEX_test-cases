
Xmin= 0.;
Xmax= 4.;
Ymin=0.;
Ymax=3.;

nbVal=40;

Xl=linspace(Xmin,Xmax,nbVal);
Yl=linspace(Ymin,Ymax,nbVal);

[Xm,Ym]=meshgrid(Xl,Yl);
vF=zeros(size(Xm));
gX=vF;
gY=vF;

for it=1:numel(Xm)
    [vF(it),gvF]=funSILEX_ZRfixed(Xm(it),Ym(it));
    gX(it)=gvF(1);
    gY(it)=gvF(2);
end

save(['paraXY_' num2str(nbVal^2) '.mat']);
