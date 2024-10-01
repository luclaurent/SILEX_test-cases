fileOPti='2016-11-25_21-04-37_optiEGO.mat';
FF=load(fileOPti);
fileOPti='2016-11-25_18-54-52_optiEGO.mat';
FFF=load(fileOPti);
close all

itO=2;
f=figure;ii=1;
for itO=1:numel(FF.optiO);
    
nbIT=size(FF.optiO{itO}.valCritEGO,1);
if FF.optiO{itO}.samplingNs==18
    ii
    
    stairs(1:10,FFF.optiO.valCritEGO(:,2),'o-','Linewidth',2);
    hold on
    stairs(1:nbIT,FF.optiO{itO}.valCritEGO(:,2),'o-','Linewidth',2);
    %ll{ii}=num2str(FF.optiO{itO}.samplingNs);
    ii=ii+1;
end
end
h=legend('$n_s=15$','$n_s=18$');
set(h,'Interpreter','latex')
line([1 10],[FF.optiO{itO}.optiTolEI FF.optiO{itO}.optiTolEI],'LineWidth',2,'Color','r');
typeYaxis='log';
set(gca, 'YScale', typeYaxis);
xlabel('Iterations EGO')
ylabel('Critère de Huang')
saveas(f,'export/EGO_crit.pdf')
saveas(f,'export/EGO_crit.eps')
matlab2tikz('standalone',true,'export/EGO_crit.tex');

f=figure;
for itO=1:numel(FF.optiO);
    
nbIT=size(FF.optiO{itO}.valCritEGO,1);
if FF.optiO{itO}.samplingNs==18
    ii
    
    stairs(1:10,FFF.optiO.minM,'o-','Linewidth',2);
    hold on
    stairs(1:nbIT,FF.optiO{itO}.minM,'o-','Linewidth',2);
end
end
h=legend('$n_s=15$','$n_s=18$');
set(h,'Interpreter','latex')
set(gca, 'YScale', typeYaxis);
xlabel('Iterations EGO')
ylabel('Minimum')

saveas(f,'export/EGO_fMin.pdf')
saveas(f,'export/EGO_fMin.eps')
matlab2tikz('standalone',true,'export/EGO_fMin.tex');