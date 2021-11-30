function ComponentIntegralErrors_Alt
% Authors: Abhishek Chatterjee, Guillaume James, and Bernard Brogliato
% Address: Univ. Grenoble Alpes, INRIA, CNRS, Grenoble INP, LJK, Grenoble
%          38000 France 
alp=3/2; 
prefixA='Figures/matfig/Components/all/'; prefixB='Figures/pdf/Components/all/';

bet=1:0.01:4;
meanAbsErr=zeros(12,length(bet));
rmsAbsErr=zeros(12,length(bet));
maxAbsErr=zeros(12,length(bet));
maxRelErr=zeros(12,length(bet));
for i=1:length(bet)
    the=0:0.01:1;
    Mstore=zeros(5,length(the)); AbsErrMstore=zeros(4,length(the)); RelErrMstore=zeros(4,length(the));
    Qstore=zeros(5,length(the)); AbsErrQstore=zeros(4,length(the)); RelErrQstore=zeros(4,length(the));
    
    Jstore=zeros(3,length(the)); AbsErrJstore=zeros(2,length(the)); RelErrJstore=zeros(2,length(the)); 
    Istore=zeros(3,length(the)); AbsErrIstore=zeros(2,length(the)); RelErrIstore=zeros(2,length(the));
    for j=1:length(the)
        [J,M]=JMnumerical(alp,bet(i),the(j));
        Jstore(1,j)=J; 
        Jstore(2,j)=Japprox(alp,bet(i),the(j),'M');
        Jstore(3,j)=Japprox(alp,bet(i),the(j),'cubic');
        
        Mstore(1,j)=M;
        Mstore(2,j)=Mapprox(alp,bet(i),the(j),'all');
        Mstore(3,j)=Mapprox(alp,bet(i),the(j),'frac');
        Mstore(4,j)=Mapprox(alp,bet(i),the(j),'xlnx');
        Mstore(5,j)=Mapprox(alp,bet(i),the(j),'x2lnx');
        
        gtil=gtilde_theta(alp,the(j));
        if the(j)==0||~isfinite(gtil)
            In=nan; Qn=nan;
        else
            [In,Qn]=VarCORCoeffsNumeric(alp,bet(i),gtil);
        end
        [Ia,Qa,Qfrac,Qxlnx,Qx2lnx,Icubic]=VarCORCoeffsApprox(alp,bet(i),the(j));
        Istore(1,j)=In; Istore(2,j)=Ia; Istore(3,j)=Icubic;
        Qstore(1,j)=Qn; Qstore(2,j)=Qa; Qstore(3,j)=Qfrac; Qstore(4,j)=Qxlnx; Qstore(5,j)=Qx2lnx;
        
        
        AbsErrJstore(1,j)=abs(Jstore(2,j) - Jstore(1,j));
        AbsErrJstore(2,j)=abs(Jstore(3,j) - Jstore(1,j));
        
        AbsErrIstore(1,j)=abs(Istore(2,j) - Istore(1,j));
        AbsErrIstore(2,j)=abs(Istore(3,j) - Istore(1,j));
       
        AbsErrMstore(1,j)=abs(Mstore(2,j) - Mstore(1,j)); 
        AbsErrMstore(2,j)=abs(Mstore(3,j) - Mstore(1,j));
        AbsErrMstore(3,j)=abs(Mstore(4,j) - Mstore(1,j));
        AbsErrMstore(4,j)=abs(Mstore(5,j) - Mstore(1,j));
        
        AbsErrQstore(1,j)=abs(Qstore(2,j) - Qstore(1,j));
        AbsErrQstore(2,j)=abs(Qstore(3,j) - Qstore(1,j));
        AbsErrQstore(3,j)=abs(Qstore(4,j) - Qstore(1,j));
        AbsErrQstore(4,j)=abs(Qstore(5,j) - Qstore(1,j));
        
        RelErrJstore(1,j)=abs((Jstore(2,j)-Jstore(1,j))/Jstore(1,j));
        RelErrJstore(2,j)=abs((Jstore(3,j)-Jstore(1,j))/Jstore(1,j));
        RelErrIstore(1,j)=abs((Istore(2,j)-Istore(1,j))/Istore(1,j));
        RelErrIstore(2,j)=abs((Istore(3,j)-Istore(1,j))/Istore(1,j));
        
        RelErrMstore(1,j)=abs((Mstore(2,j)-Mstore(1,j))/Mstore(1,j));
        RelErrMstore(2,j)=abs((Mstore(3,j)-Mstore(1,j))/Mstore(1,j));
        RelErrMstore(3,j)=abs((Mstore(4,j)-Mstore(1,j))/Mstore(1,j));
        RelErrMstore(4,j)=abs((Mstore(5,j)-Mstore(1,j))/Mstore(1,j));
        
        RelErrQstore(1,j)=abs((Qstore(2,j)-Qstore(1,j))/Qstore(1,j));
        RelErrQstore(2,j)=abs((Qstore(3,j)-Qstore(1,j))/Qstore(1,j));
        RelErrQstore(3,j)=abs((Qstore(4,j)-Qstore(1,j))/Qstore(1,j));
        RelErrQstore(4,j)=abs((Qstore(5,j)-Qstore(1,j))/Qstore(1,j));
        
        
    end 
    
    % Store Maximum Relative Errors
    maxRelErr(1,i)=max(RelErrJstore(1,:));
    maxRelErr(2,i)=max(RelErrJstore(2,:));
    
    maxRelErr(3,i)=max(RelErrIstore(1,:));
    maxRelErr(4,i)=max(RelErrIstore(2,:));
    
    maxRelErr(5,i)=max(RelErrMstore(1,:));
    maxRelErr(6,i)=max(RelErrMstore(2,:));
    maxRelErr(7,i)=max(RelErrMstore(3,:));
    maxRelErr(8,i)=max(RelErrMstore(4,:));
    
    maxRelErr(9,i)=max(RelErrQstore(1,:));
    maxRelErr(10,i)=max(RelErrQstore(2,:));
    maxRelErr(11,i)=max(RelErrQstore(3,:));
    maxRelErr(12,i)=max(RelErrQstore(4,:));
    
    % Store Mean Absolute Errors
    meanAbsErr(1,i)=mean(AbsErrJstore(1,:));
    meanAbsErr(2,i)=mean(AbsErrJstore(2,:));
    
    meanAbsErr(3,i)=mean(AbsErrIstore(1,:));
    meanAbsErr(4,i)=mean(AbsErrIstore(2,:));
    
    meanAbsErr(5,i)=mean(AbsErrMstore(1,:));
    meanAbsErr(6,i)=mean(AbsErrMstore(2,:));
    meanAbsErr(7,i)=mean(AbsErrMstore(3,:));
    meanAbsErr(8,i)=mean(AbsErrMstore(4,:));
    
    meanAbsErr(9,i)=mean(AbsErrQstore(1,:));
    meanAbsErr(10,i)=mean(AbsErrQstore(2,:));
    meanAbsErr(11,i)=mean(AbsErrQstore(3,:));
    meanAbsErr(12,i)=mean(AbsErrQstore(4,:));
    
    % Store RMS of Absolute Errors 
    rmsAbsErr(1,i)=sqrt(mean(AbsErrJstore(1,:).^2));
    rmsAbsErr(2,i)=sqrt(mean(AbsErrJstore(2,:).^2));
    
    rmsAbsErr(3,i)=sqrt(mean(AbsErrIstore(1,:).^2));
    rmsAbsErr(4,i)=sqrt(mean(AbsErrIstore(2,:).^2));
    
    rmsAbsErr(5,i)=sqrt(mean(AbsErrMstore(1,:).^2));
    rmsAbsErr(6,i)=sqrt(mean(AbsErrMstore(2,:).^2));
    rmsAbsErr(7,i)=sqrt(mean(AbsErrMstore(3,:).^2));
    rmsAbsErr(8,i)=sqrt(mean(AbsErrMstore(4,:).^2));
    
    rmsAbsErr(9,i)=sqrt(mean(AbsErrQstore(1,:).^2));
    rmsAbsErr(10,i)=sqrt(mean(AbsErrQstore(2,:).^2));
    rmsAbsErr(11,i)=sqrt(mean(AbsErrQstore(3,:).^2));
    rmsAbsErr(12,i)=sqrt(mean(AbsErrQstore(4,:).^2));
    
    % Store Max of Absolute Errors 
    maxAbsErr(1,i)=max(AbsErrJstore(1,:));
    maxAbsErr(2,i)=max(AbsErrJstore(2,:));
    
    maxAbsErr(3,i)=max(AbsErrIstore(1,:));
    maxAbsErr(4,i)=max(AbsErrIstore(2,:));
    
    maxAbsErr(5,i)=max(AbsErrMstore(1,:));
    maxAbsErr(6,i)=max(AbsErrMstore(2,:));
    maxAbsErr(7,i)=max(AbsErrMstore(3,:));
    maxAbsErr(8,i)=max(AbsErrMstore(4,:));
    
    maxAbsErr(9,i)=max(AbsErrQstore(1,:));
    maxAbsErr(10,i)=max(AbsErrQstore(2,:));
    maxAbsErr(11,i)=max(AbsErrQstore(3,:));
    maxAbsErr(12,i)=max(AbsErrQstore(4,:));
end

% Plots for Maximum Relative Error
figure('Color','w')
semilogy(bet,maxRelErr(1,:),'-b','LineWidth',2); 
%legend('$M_{\alpha,\beta} (\theta) = a_0 \theta^{\beta-1/2} + a_1 \theta \ln \theta  + a_2 \theta^2 \ln \theta + a_3 \theta^2 + a_4 \theta + a_5$',,'Interpreter','latex','FontSize',16)
title('Max. Relative Error of $J_{\alpha,\beta}(\theta)$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Max. Relative Error, $J_{\alpha,\beta}(\theta)$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'Alt_JmaxRelErr.fig'));
print(strcat(prefixB,'Alt_JmaxRelErr.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,maxRelErr(1,:),'-b','LineWidth',2); 
hold on 
semilogx(bet,maxRelErr(2,:),'--g','LineWidth',2);
legend('J_{\alpha,\beta} (\theta) = \frac{\alpha+1}{2\beta} M_{\alpha,\alpha + \beta+1} (\theta) + \frac{\theta - 1}{2\beta} M_{\alpha,\beta+1} (\theta)',...
    'J_{\alpha,\beta} (\theta) = a \theta^3 + b \theta^2 + c \theta + d','Interpreter','latex','FontSize',16)
title('Max. Relative Error of $J_{\alpha,\beta}(\theta)$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Max. Relative Error, $J_{\alpha,\beta}(\theta)$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'Alt_CompJmaxRelErr.fig'));
print(strcat(prefixB,'Alt_CompJmaxRelErr.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,maxRelErr(3,:),'-b','LineWidth',2); 
title('Max. Relative Error of $\mathcal{I}(\tilde{g}(\theta))$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Max. Relative Error, $\mathcal{I}(\tilde{g}(\theta))$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'Alt_ImaxRelErr.fig'));
print(strcat(prefixB,'Alt_ImaxRelErr.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,maxRelErr(3,:),'-b','LineWidth',2); 
hold on
semilogy(bet,maxRelErr(4,:),'--g','LineWidth',2);
legend('Rewrite','Cubic')
title('Max. Relative Error of $\mathcal{I}(\tilde{g}(\theta))$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Max. Relative Error, $\mathcal{I}(\tilde{g}(\theta))$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'Alt_CompImaxRelErr.fig'));
print(strcat(prefixB,'Alt_CompImaxRelErr.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,maxRelErr(5,:),'-b','LineWidth',2); 
title('Max. Relative Error of $M_{\alpha,\beta}(\theta)$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Max. Relative Error, $M_{\alpha,\beta}(\theta)$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'Alt_MmaxRelErr.fig'));
print(strcat(prefixB,'Alt_MmaxRelErr.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,maxRelErr(5,:),'-b','LineWidth',2); 
hold on
semilogy(bet,maxRelErr(6,:),'-r');
hold on 
semilogy(bet,maxRelErr(7,:),'-g');
hold on
semilogy(bet,maxRelErr(8,:),'-m');
legend('$M_{\alpha,\beta} (\theta)$','$M^q_{\alpha,\beta} (\theta)$','$M^r_{\alpha,\beta} (\theta)$', '$M^s_{\alpha,\beta} (\theta)$', 'Interpreter', 'latex')
title('Compare Max. Relative Error of $M_{\alpha,\beta}(\theta)$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Compare Max. Relative Error, $M_{\alpha,\beta}(\theta)$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'Alt_CompMmaxRelErr.fig'));
print(strcat(prefixB,'Alt_CompMmaxRelErr.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,maxRelErr(9,:),'-b','LineWidth',2); 
title('Max. Relative Error of $\mathcal{Q}(\tilde{g}(\theta))$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Max. Relative Error, $\mathcal{Q}(\tilde{g}(\theta))$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'Alt_QmaxRelErr.fig'));
print(strcat(prefixB,'Alt_QmaxRelErr.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,maxRelErr(9,:),'-b','LineWidth',2);
hold on
semilogy(bet,maxRelErr(10,:),'--r');
hold on 
semilogy(bet,maxRelErr(11,:),'--g');
hold on
semilogy(bet,maxRelErr(12,:),'--m');
legend('$Piecewise$','$\theta^{\beta-1/2}$','$\theta \ln{\theta}$', '$\theta^2 \ln \theta$', 'Interpreter', 'latex')
title('Compare Max. Relative Error of $\mathcal{Q}(\tilde{g}(\theta))$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Max. Relative Error, $\mathcal{Q}(\tilde{g}(\theta))$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'Alt_CompQmaxRelErr.fig'));
print(strcat(prefixB,'Alt_CompQmaxRelErr.pdf'),'-dpdf');

% Plots of Mean Absolute Error
%{
figure('Color','w')
semilogy(bet,meanAbsErr(1,:),'-b','LineWidth',2); 
title('Mean of the Absolute Errors, $J_{\alpha,\beta}(\theta)$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Mean of the Absolute Errors, $J_{\alpha,\beta}(\theta)$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'JmeanAbsErr.fig'));
print(strcat(prefixB,'JmeanAbs.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,meanAbsErr(1,:),'-b','LineWidth',2);
title('Mean of the Absolute Errors, $J_{\alpha,\beta}(\theta)$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Mean of the Absolute Errors, $J_{\alpha,\beta}(\theta)$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'JmeanAbsErr.fig'));
print(strcat(prefixB,'JmeanAbs.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,meanAbsErr(2,:),'-b','LineWidth',2); 
title('Mean of the Absolute Errors, $\mathcal{I}(\tilde{g}(\theta))$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Mean of the Absolute Errors, $\mathcal{I}(\tilde{g}(\theta))$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'ImeanAbsErr.fig'));
print(strcat(prefixB,'ImeanAbsErr.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,meanAbsErr(3,:),'-b','LineWidth',2); 
title('Mean of the Absolute Errors, $M_{\alpha,\beta}(\theta)$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Mean of the Absolute Errors, $M_{\alpha,\beta}(\theta)$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'MmeanAbsErr.fig'));
print(strcat(prefixB,'MmeanAbsErr.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,meanAbsErr(3,:),'-b','LineWidth',2); 
hold on
semilogy(bet,meanAbsErr(4,:),'--r');
hold on 
semilogy(bet,meanAbsErr(5,:),'--g');
hold on
semilogy(bet,meanAbsErr(6,:),'--m');
legend('$Piecewise$','$\theta^{\beta-1/2}$','$\theta \ln{\theta}$', '$\theta^2 \ln \theta$', 'Interpreter', 'latex')
title('Compare Mean of the Absolute Errors of $M_{\alpha,\beta}(\theta)$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Mean of the Absolute Errors, $M_{\alpha,\beta}(\theta)$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'CompMmeanAbsErr.fig'));
print(strcat(prefixB,'CompMmeanAbsErr.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,meanAbsErr(7,:),'-b','LineWidth',2); 
title('Mean of the Absolute Errors, $\mathcal{Q}(\tilde{g}(\theta))$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Mean of the Absolute Errors, $\mathcal{Q}(\tilde{g}(\theta))$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'QmeanAbsErr.fig'));
print(strcat(prefixB,'QmeanAbsErr.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,meanAbsErr(7,:),'-b','LineWidth',2);
hold on
semilogy(bet,meanAbsErr(8,:),'--r');
hold on 
semilogy(bet,meanAbsErr(9,:),'--g');
hold on
semilogy(bet,meanAbsErr(10,:),'--m');
legend('$Piecewise$','$\theta^{\beta-1/2}$','$\theta \ln{\theta}$', '$\theta^2 \ln \theta$', 'Interpreter', 'latex')
title('Compare Mean of the Absolute Errors of $\mathcal{Q}(\tilde{g}(\theta))$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Mean of the Absolute Errors, $\mathcal{Q}(\tilde{g}(\theta))$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'CompQmeanAbsErr.fig'));
print(strcat(prefixB,'CompQmeanAbsErr.pdf'),'-dpdf');


% Plots of Root Mean Square of the Absolute Error

figure('Color','w')
semilogy(bet,rmsAbsErr(1,:),'-b','LineWidth',2); 
title('RMS of the Absolute Errors, $J_{\alpha,\beta}(\theta)$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('RMS of the Absolute Errors, $J_{\alpha,\beta}(\theta)$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'JrmsAbsErr.fig'));
print(strcat(prefixB,'JrmsAbs.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,rmsAbsErr(2,:),'-b','LineWidth',2); 
title('RMS of the Absolute Errors, $\mathcal{I}(\tilde{g}(\theta))$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('RMS of the Absolute Errors, $\mathcal{I}(\tilde{g}(\theta))$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'IrmsAbsErr.fig'));
print(strcat(prefixB,'IrmsAbsErr.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,rmsAbsErr(3,:),'-b','LineWidth',2); 
title('RMS of the Absolute Errors, $M_{\alpha,\beta}(\theta)$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('RMS of the Absolute Errors, $M_{\alpha,\beta}(\theta)$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'MrmsAbsErr.fig'));
print(strcat(prefixB,'MrmsAbsErr.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,rmsAbsErr(3,:),'-b','LineWidth',2); 
hold on
semilogy(bet,rmsAbsErr(4,:),'--r');
hold on 
semilogy(bet,rmsAbsErr(5,:),'--g');
hold on
semilogy(bet,rmsAbsErr(6,:),'--m');
legend('$Piecewise$','$\theta^{\beta-1/2}$','$\theta \ln{\theta}$', '$\theta^2 \ln \theta$', 'Interpreter', 'latex')
title('Compare RMS of the Absolute Errors of $M_{\alpha,\beta}(\theta)$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('RMS of the Absolute Errors, $M_{\alpha,\beta}(\theta)$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'CompMrmsAbsErr.fig'));
print(strcat(prefixB,'CompMrmsAbsErr.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,rmsAbsErr(7,:),'-b','LineWidth',2); 
title('RMS of the Absolute Errors, $\mathcal{Q}(\tilde{g}(\theta))$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('RMS of the Absolute Errors, $\mathcal{Q}(\tilde{g}(\theta))$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'QrmsAbsErr.fig'));
print(strcat(prefixB,'QrmsAbsErr.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,rmsAbsErr(7,:),'-b','LineWidth',2);
hold on
semilogy(bet,rmsAbsErr(8,:),'--r');
hold on 
semilogy(bet,rmsAbsErr(9,:),'--g');
hold on
semilogy(bet,rmsAbsErr(10,:),'--m');
legend('$Piecewise$','$\theta^{\beta-1/2}$','$\theta \ln{\theta}$', '$\theta^2 \ln \theta$', 'Interpreter', 'latex')
title('Compare RMS of the Absolute Errors of $\mathcal{Q}(\tilde{g}(\theta))$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('RMS of the Absolute Errors, $\mathcal{Q}(\tilde{g}(\theta))$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'CompQrmsAbsErr.fig'));
print(strcat(prefixB,'CompQrmsAbsErr.pdf'),'-dpdf');

% Plots of Max. Absolute Error

figure('Color','w')
semilogy(bet,maxAbsErr(1,:),'-b','LineWidth',2); 
title('Max. of the Absolute Errors, $J_{\alpha,\beta}(\theta)$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Max. of the Absolute Errors, $J_{\alpha,\beta}(\theta)$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'JmaxAbsErr.fig'));
print(strcat(prefixB,'JmaxAbs.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,maxAbsErr(2,:),'-b','LineWidth',2); 
title('Max. of the Absolute Errors, $\mathcal{I}(\tilde{g}(\theta))$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Max. of the Absolute Errors, $\mathcal{I}(\tilde{g}(\theta))$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'ImaxAbsErr.fig'));
print(strcat(prefixB,'ImaxAbsErr.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,maxAbsErr(3,:),'-b','LineWidth',2); 
title('Max. of the Absolute Errors, $M_{\alpha,\beta}(\theta)$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Max. of the Absolute Errors, $M_{\alpha,\beta}(\theta)$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'MmaxAbsErr.fig'));
print(strcat(prefixB,'MmaxAbsErr.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,maxAbsErr(3,:),'-b','LineWidth',2); 
hold on
semilogy(bet,maxAbsErr(4,:),'--r');
hold on 
semilogy(bet,maxAbsErr(5,:),'--g');
hold on
semilogy(bet,maxAbsErr(6,:),'--m');
legend('$Piecewise$','$\theta^{\beta-1/2}$','$\theta \ln{\theta}$', '$\theta^2 \ln \theta$', 'Interpreter', 'latex')
title('Compare Max. of the Absolute Errors of $M_{\alpha,\beta}(\theta)$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Max. of the Absolute Errors, $M_{\alpha,\beta}(\theta)$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'CompMmaxAbsErr.fig'));
print(strcat(prefixB,'CompMmaxAbsErr.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,maxAbsErr(7,:),'-b','LineWidth',2); 
title('Max. of the Absolute Errors, $\mathcal{Q}(\tilde{g}(\theta))$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Max. of the Absolute Errors, $\mathcal{Q}(\tilde{g}(\theta))$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'QmaxAbsErr.fig'));
print(strcat(prefixB,'QmaxAbsErr.pdf'),'-dpdf');

figure('Color','w')
semilogy(bet,maxAbsErr(7,:),'-b','LineWidth',2);
hold on
semilogy(bet,maxAbsErr(8,:),'--r');
hold on 
semilogy(bet,maxAbsErr(9,:),'--g');
hold on
semilogy(bet,maxAbsErr(10,:),'--m');
legend('$Piecewise$','$\theta^{\beta-1/2}$','$\theta \ln{\theta}$', '$\theta^2 \ln \theta$', 'Interpreter', 'latex')
title('Compare Max. of the Absolute Errors of $\mathcal{Q}(\tilde{g}(\theta))$ w.r.t. $\beta$','Interpreter','latex','FontSize',16)
ylabel('Max. of the Absolute Errors, $\mathcal{Q}(\tilde{g}(\theta))$','Interpreter','latex','FontSize',16)
xlabel('$\beta$','Interpreter','latex','FontSize',16)
savefig(strcat(prefixA,'CompQmaxAbsErr.fig'));
print(strcat(prefixB,'CompQmaxAbsErr.pdf'),'-dpdf');

%}
end