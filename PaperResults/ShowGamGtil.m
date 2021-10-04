function ShowGamGtil
%v0=0:0.0001:0.1;
v0=0:0.0001:1;
alp=3/2; bet=3/2;
%gam0=0.001;m=2.05e-3; k=9.858e9; 
[k,m,gam0] = get_params;
gam=gam0.*(v0.^((2*bet)/(alp+1) - 1)).*((k/m)^(1-(bet/(alp+1))));
figure('color','w')
plot(v0,gam,'-b','LineWidth',1);
title('$\gamma$ vs. $v_0$, given $\alpha=\beta=3/2$ and $\gamma_0=1.5237 \times 10^{-6}$','Interpreter','Latex','FontSize',16)
xlabel('$v_0$','Interpreter','Latex','FontSize',16)
ylabel('$\gamma$','Interpreter','Latex','FontSize',16)
savefig('Figures/matfig/ParamView/gammaWRTv0.fig')
print('Figures/pdf/ParamView/gammaWRTv0.pdf','-dpdf')

g=9.8;
%f=-0.1:0.0001:0.1;
f=-2:0.0001:2;
v=logspace(-1,0,5);
ls={'-b','-r','--g','-.m',':k'}; leg={};
figure('Color','w')
for i=1:length(v)
    gtil=((m/k)^(1/(alp+1)))*(v(i)^(-2*alp/(alp+1))).*(g+f./m);
    plot(f,gtil,ls{i},'LineWidth',1);
    leg{i}=strcat('$v_0 =',num2str(v(i),'%.3f'),' {\rm m/s} $');
    hold on
end
title('$\tilde{g}$ vs. $F$, given $\alpha=\beta=3/2$ and $\gamma_0=1.5237 \times 10^{-6}$','Interpreter','Latex','FontSize',16)
xlabel('$F$','Interpreter','Latex','FontSize',16)
ylabel('$\tilde{g}$','Interpreter','Latex','FontSize',16)
legend(leg{:},'Interpreter','Latex','FontSize',12)
ax=gca; 
ax.YAxis.Exponent=0;
savefig('Figures/matfig/ParamView/gtilWRTF.fig')
print('Figures/pdf/ParamView/gtilWRTF.pdf','-dpdf')
end

function [k,m,gam0]=get_params
e=0.893; %steel; (from Kuwabara-Kono paper)
gam = 0.438*(1-e); 
v0 = 0.5; 
E=21.1e10; R=1.65e-2; m=154e-3; 
k=(4/3)*sqrt(R)*E;
gam0 = gam*(v0^(-1/5))*((k/m)^(-2/5));
end