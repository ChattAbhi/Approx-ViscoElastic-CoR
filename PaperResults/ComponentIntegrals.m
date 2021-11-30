function ComponentIntegrals(choice)
% Authors: Abhishek Chatterjee, Guillaume James, and Bernard Brogliato
% Address: Univ. Grenoble Alpes, INRIA, CNRS, Grenoble INP, LJK, Grenoble
%          38000 France 
alp=3/2; 

if choice==1
    bet=1; prefixA='Figures/matfig/Components/case_1/'; prefixB='Figures/pdf/Components/case_1/';
elseif choice==2
    bet=5/4; prefixA='Figures/matfig/Components/case_2/'; prefixB='Figures/pdf/Components/case_2/';
elseif choice==3
    bet=3/2; prefixA='Figures/matfig/Components/case_3/'; prefixB='Figures/pdf/Components/case_3/';
elseif choice==4
    bet=5/2; prefixA='Figures/matfig/Components/case_4/'; prefixB='Figures/pdf/Components/case_4/';
else
    bet=2.51; prefixA='Figures/matfig/Components/tests'; prefixB='Figures/pdf/Components/tests';
end

the=0:0.01:1;
Jstore=zeros(2,length(the)); ErJ=zeros(1,length(the));
Mstore=zeros(2,length(the)); ErM=zeros(1,length(the));
Istore=zeros(2,length(the)); ErI=zeros(1,length(the));
Qstore=zeros(2,length(the)); ErQ=zeros(1,length(the));
gtilde=zeros(1,length(the));
for i=1:length(the)
    [J,M]=JMnumerical(alp,bet,the(i));
    Jstore(1,i)=J; Mstore(1,i)=M;
    Jstore(2,i)=Japprox(alp,bet,the(i));
    Mstore(2,i)=Mapprox(alp,bet,the(i));
    ErJ(1,i)=abs((Jstore(2,i)-Jstore(1,i))/Jstore(1,i));
    ErM(i)=abs((Mstore(2,i)-Mstore(1,i))/Mstore(1,i));
    
    gtilde(i)=gtilde_theta(alp,the(i));
    if the(i)==0||~isfinite(gtilde(i))
        I1=nan; Q1=nan;
    else
        [I1,Q1]=VarCORCoeffsNumeric(alp,bet,gtilde(i));
    end
    [I2,Q2]=VarCORCoeffsApprox(alp,bet,the(i));
    Istore(1,i)=I1; Qstore(1,i)=Q1;
    Istore(2,i)=I2; Qstore(2,i)=Q2;
    ErI(i)=abs((Istore(2,i)-Istore(1,i))/Istore(1,i));
    ErQ(i)=abs((Qstore(2,i)-Qstore(1,i))/Qstore(1,i));
end

figure('Color','w')
plot(the,Jstore(1,:),'-b','LineWidth',1)
hold on
plot(the,Jstore(2,:),'--r','LineWidth',1)
hold off 
title(strcat('$J_{\alpha,\beta}(\theta)$, with $\alpha = ',num2str(alp),'$ and $\beta =',num2str(bet),'$'),'Interpreter','latex','FontSize',16)
legend('Numerical','Approximation','Interpreter','latex','FontSize',16,'Location','Best')
ylabel('$J_{\alpha,\beta}(\theta)$','Interpreter','latex','FontSize',20)
xlabel('$\theta$','Interpreter','latex','FontSize',20)
savefig(strcat(prefixA,'Japprox.fig'));
print(strcat(prefixB,'Japprox.pdf'),'-dpdf');


figure('Color','w')
plot(the,ErJ,'--r')
title(strcat('Relative Error of $J_{\alpha,\beta}(\theta)$, with $\alpha =',num2str(alp),'$ and $\beta =',num2str(bet),'$'),'Interpreter','latex','FontSize',16)
%legend('Approximation','Interpreter','latex','FontSize',16,'Location','Best')
ylabel('Relative Error','Interpreter','latex','FontSize',20)
xlabel('$\theta$','Interpreter','latex','FontSize',20)
savefig(strcat(prefixA,'JapproxErr.fig'));
print(strcat(prefixB,'JapproxErr.pdf'),'-dpdf');

figure('Color','w')
plot(the,Mstore(1,:),'-b','LineWidth',1)
hold on
plot(the,Mstore(2,:),'--r','LineWidth',1)
hold off 
title(strcat('$M_{\alpha,\beta}(\theta)$, given $\alpha = ',num2str(alp),'$ and $\beta =',num2str(bet),'$'),'Interpreter','latex','FontSize',16)
legend('Numerical','Analytical','Interpreter','latex','FontSize',16,'Location','Best')
ylabel('$M_{\alpha,\beta}(\theta)$','Interpreter','latex','FontSize',20)
xlabel('$\theta$','Interpreter','latex','FontSize',20)
savefig(strcat(prefixA,'MCapprox.fig'));
print(strcat(prefixB,'MCapprox.pdf'),'-dpdf');

figure('Color','w')
plot(the,ErM,'-b')
title(strcat('Relative Error of $M_{\alpha,\beta}(\theta)$, with $\alpha = ', num2str(alp),'$ and $\beta =',num2str(bet),'$'),'Interpreter','latex','FontSize',16)
ylabel('Relative Error','Interpreter','latex','FontSize',20)
xlabel('$\theta$','Interpreter','latex','FontSize',20)
savefig(strcat(prefixA,'MCapproxErr.fig'));
print(strcat(prefixB,'MCapproxErr.pdf'),'-dpdf');

figure('Color','w')
plot(gtilde,Istore(1,:),'-b','LineWidth',1)
hold on 
plot(gtilde,Istore(2,:),'--r','LineWidth',1)
title(strcat('$\mathcal{I}(\tilde{g}(\theta))=K(\theta)J(\theta)$, with $\alpha =', num2str(alp),'$ and $\beta =',num2str(bet),'$'),'Interpreter','latex','FontSize',16)
legend('Numerical','Analytical','Interpreter','latex','FontSize',16,'Location','Best')
ylabel('$\mathcal{I}(\tilde{g}(\theta))$','Interpreter','latex','FontSize',20)
xlabel('$\tilde{g}(\theta)$','Interpreter','latex','FontSize',20)
savefig(strcat(prefixA,'Iapprox.fig'));
print(strcat(prefixB,'Iapprox.pdf'),'-dpdf');

figure('Color','w')
plot(gtilde,ErI,'--r')
title(strcat('Relative Error of $\mathcal{I}(\tilde{g}(\theta))=K(\theta)J(\theta)$, with $\alpha = ', num2str(alp),'$ and $\beta =',num2str(bet),'$'),'Interpreter','latex','FontSize',16)
ylabel('Relative Error','Interpreter','latex','FontSize',20)
xlabel('$\tilde{g}(\theta)$','Interpreter','latex','FontSize',20)
savefig(strcat(prefixA,'IapproxErr.fig'));
print(strcat(prefixB,'IapproxErr.pdf'),'-dpdf');

figure('Color','w')
plot(gtilde,Qstore(1,:),'-b','LineWidth',1)
hold on 
plot(gtilde,Qstore(2,:),'--r','LineWidth',1)
title(strcat('$\mathcal{Q}(\tilde{g}(\theta))=L(\theta)M_{\alpha,\beta}(\theta)$, with $\alpha =', num2str(alp),'$ and $\beta =',num2str(bet),'$'),'Interpreter','latex','FontSize',16)
legend('Numerical','Analytical','Interpreter','latex','FontSize',16,'Location','Best')
ylabel('$\mathcal{Q}(\tilde{g}(\theta))$','Interpreter','latex','FontSize',20)
xlabel('$\tilde{g}(\theta)$','Interpreter','latex','FontSize',20)
savefig(strcat(prefixA,'QCapprox.fig'));
print(strcat(prefixB,'QCapprox.pdf'),'-dpdf');

figure('Color','w')
plot(gtilde,ErQ,'--r')
title(strcat('Relative Error of $\mathcal{Q}(\tilde{g}(\theta))=L(\theta) M_{\alpha,\beta}(\theta)$, with $\alpha = ', num2str(alp),'$ and $\beta =',num2str(bet),'$'),'Interpreter','latex','FontSize',16)
ylabel('Relative Error','Interpreter','latex','FontSize',20)
xlabel('$\tilde{g}(\theta)$','Interpreter','latex','FontSize',20)
savefig(strcat(prefixA,'QCapproxErr.fig'));
print(strcat(prefixB,'QCapproxErr.pdf'),'-dpdf');

end