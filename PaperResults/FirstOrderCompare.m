function FirstOrderCompare(gam)
alp =1.5; %alpha
bet =1.5; %beta

m=0.00205; k=6971600000; 
v0=0.246; g=9.8;
gtil_val=((m/k)^(1/(alp+1)))*(v0^((-2*alp)/(alp+1)))*g;

gtilde_min = 0; 
gtilde_max = 10; 
num_steps=100;

gtilde_range = linspace(gtilde_min,gtilde_max,num_steps);
E_num=zeros(1,length(gtilde_range));
E_anl=zeros(1,length(gtilde_range));
[C0,C1,~]=ConstCORCoeffs(alp,bet);
for i=1:length(gtilde_range)
    E_num(i)=numericCOR(alp,bet,gam,gtilde_range(i));
    E_anl(i)=1 - gam*C0 - gam*gtilde_range(i)*C1;
end

figure('Color','w')
plot(gtilde_range,E_num,'-b','LineWidth',1)
hold on
plot(gtilde_range,E_anl,'--r','LineWidth',1)
legend('Num. Integ., [Eq.~(4)]','Approx. (Lin. $\tilde{g}$)  - $O(\gamma)$, [Eq.~(22)]','Interpreter','latex','FontSize',16,'Location','Best')
hold off
title(strcat('$e$ w.r.t change in $\tilde{g}$ ($\gamma$ =',num2str(gam),')'),'Interpreter','latex','FontSize',20)
xlabel('$$\tilde{g}$$','Interpreter','latex','FontSize',20)
ylabel('$$e$$','Interpreter','latex','FontSize',20)

gtilde_min = 0; 
gtilde_max = 1; 
num_steps=8;

gtilde_range = linspace(gtilde_min,gtilde_max,num_steps);
E_num=zeros(1,length(gtilde_range));
E_anl=zeros(1,length(gtilde_range));
for i=1:length(gtilde_range)
    E_num(i)=numericCOR(alp,bet,gam,gtilde_range(i));
    E_anl(i) =1 - gam*C0 - gam*gtilde_range(i)*C1;
end
EOrgRange=[E_num(1),E_num(end)];

fig=figure('Color','w');
ax_base=axes(fig);
h(1)=plot(ax_base,gtilde_range,E_num,'-b','LineWidth',1);
hold on
h(2)=plot(ax_base,gtilde_range,E_anl,'--r','LineWidth',1);
side=0.0003; top=0.00005; maxE=max([E_num,E_anl]); fact1=100; fact2=10;
small_box=[0,(1+top)*maxE;fact1*side,(1+top)*maxE;fact1*side,(1-fact2*side+top)*maxE;0,(1-fact2*side+top)*maxE;0,(1+top)*maxE]';
%h(3)=plot(small_box(1,:),small_box(2,:),'-k');
legend(h,{'Num. Integ., [Eq.~(4)]','Approx. (Lin. $\tilde{g}$)  - $O(\gamma)$, [Eq.~(22)]'},'Interpreter','latex','FontSize',16,'Location','Best')
hold off
title(strcat('$e$ w.r.t change in $\tilde{g}$ ($\gamma$ =',num2str(gam),')'),'Interpreter','latex','FontSize',20)
xlabel('$$\tilde{g}$$','Interpreter','latex','FontSize',20)
ylabel('$$e$$','Interpreter','latex','FontSize',20)

gtilde_min = 0; 
gtilde_max = 4*gtil_val; 
num_steps=10;

gtilde_range = linspace(gtilde_min,gtilde_max,num_steps);
E_num=zeros(1,length(gtilde_range));
E_anl=zeros(1,length(gtilde_range));
for i=1:length(gtilde_range)
    E_num(i)=numericCOR(alp,bet,gam,gtilde_range(i));
    E_anl(i)=1 - gam*C0 - gam*gtilde_range(i)*C1;
end
pos=get(ax_base,'Position');
scl=0.65; dims=[pos(3)-pos(1),pos(4)-pos(2)];
ax_zoom=axes(fig,'OuterPosition',[pos(1),pos(2),scl*dims(1),scl*dims(2)],'Box','on');
h(1)=plot(ax_zoom,gtilde_range,E_num,'-b','LineWidth',1);
hold on
h(2)=plot(ax_zoom,gtilde_range,E_anl,'--r','LineWidth',1);
hold on

h(3)=plot(ax_zoom,[gtil_val,gtil_val],[(1+top)*maxE,(1-side+top)*maxE],'-.k');
axis(ax_zoom,[0,4*gtil_val,(1-side+top)*maxE,(1+top)*maxE])

for i=1:length(gtilde_range)
    if gtilde_range(i)>gtil_val
        Eng0 = E_num(i-1) + ((E_num(i)- E_num(i-1))/(gtilde_range(i) - gtilde_range(i-1)))*(gtil_val - gtilde_range(i));
    end
end

Eag0=1 - gam*C0 - gam*gtil_val*C1;

AbsErrEg0=abs(Eng0-Eag0)
text(ax_base,0.6,0.9*EOrgRange(1),strcat('$\tilde{g}_0 =',num2str(round(gtil_val,4,'decimal')),'$'),'Interpreter','latex','FontSize',16)
text(ax_base,0.65,1.2*EOrgRange(2),num2str(round(AbsErrEg0,4,'decimal')),'Interpreter','latex','FontSize',16)
fprintf('The Absolute Error around g_0 = %.4f, (when gamma = %.4f) is %.4f \n',gtil_val,gam,AbsErrEg0);
end