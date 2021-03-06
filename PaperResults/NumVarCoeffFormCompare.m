function NumVarCoeffFormCompare

% Authors: Abhishek Chatterjee, Guillaume James, and Bernard Brogliato
% Address: Univ. Grenoble Alpes, INRIA, CNRS, Grenoble INP, LJK, Grenoble
%          38000 France 

%% Compare O(gamma^2) between Schwager-Poschel and numerically computed coeffeicients
alp=3/2; bet=3/2;
% Get results for specific gtilde values 
gtil_choice=[0,0.05,1,10];
z_fact=[3/1000,1/10,1/10000,1/1000];
a=0.30522;
max_gam=zeros(1,length(gtil_choice));
for i=1:length(gtil_choice)
    max_gam(i)=(a/gtil_choice(i))^(11/6);
end

prefixA='Figures/matfig/NumCoeff/'; prefixB='Figures/pdf/NumCoeff/';
num_step=1000;
gam_lim=0.1; 
if bet==3/2
    [C0,C1,C2]=ConstCORCoeffs(alp,bet);
else
    C2=0; %Order 2 constant coefficient is not known when
    [C0,C1]=ConstCORCoeffs(alp,bet);
end
for i=1:length(gtil_choice)
    E=zeros(3,num_step); Ez=zeros(3,num_step);
    R=zeros(2,num_step); Rz=zeros(2,num_step);
    A=zeros(2,num_step); Az=zeros(2,num_step);
    
    gam_range=linspace(0,min([gam_lim,1.1*max_gam(i)]),num_step);
    %gam_range_z=linspace(0,z_fact(i)*min([gam_lim,1.1*max_gam(i)]),num_step);
    
    [I,Q]=VarCORCoeffsNumeric(alp,bet,gtil_choice(i));
    gam_z_max = 2/(bet*(Q-(1/2)*I));
    if gam_z_max>0 && gam_z_max<0.1*min([gam_lim,1.1*max_gam(i)])
        gam_range_z=linspace(0,2*gam_z_max,num_step);
    else
        gam_range_z=linspace(0,z_fact(i)*min([gam_lim,1.1*max_gam(i)]),num_step);
    end
    
    for j=1:length(gam_range)
        E(1,j)= numericCOR(alp,bet,gam_range(j),gtil_choice(i));
        E(2,j)= 1 - gam_range(j)*C0 - gam_range(j)*gtil_choice(i)*C1 - C2*gam_range(j)^2; %O(gamma^2) COR with Schwager-Poschel
        esqr=1 - 2*bet*I*gam_range(j) + 2*(bet^2)*I*Q*(gam_range(j)^2);
        E(3,j)=sqrt(max([esqr,0])); %O(gamma^2) COR, as \sqrt(e^2), with numerical I(gtilde) and Q(gtilde)
        E(4,j)=1 - bet*I*gam_range(j) + (bet^2)*I*(Q - (1/2)*I)*(gam_range(j)^2); %O(gamma^2) COR, as e= 1 + h/2 + h^2/8, with numerical I(gtilde) and Q(gtilde)
        
        
        A(1,j)= abs(E(2,j) - E(1,j)); %O(gamma^2) COR with Schwager-Poschel 
        A(2,j)= abs(E(3,j) - E(1,j)); %O(gamma^2) COR, as \sqrt(e^2), with numerical I(gtilde) and Q(gtilde)
        A(3,j)= abs(E(4,j) - E(1,j)); %O(gamma^2) COR, as e= 1 + h/2 + h^2/8, with numerical I(gtilde) and Q(gtilde)
        
        R(1,j)= abs((E(2,j) - E(1,j))/E(1,j)); %O(gamma^2) COR with Schwager-Poschel 
        R(2,j)= abs((E(3,j) - E(1,j))/E(1,j)); %O(gamma^2) COR, as \sqrt(e^2), with numerical I(gtilde) and Q(gtilde)
        R(3,j)= abs((E(4,j) - E(1,j))/E(1,j)); %O(gamma^2) COR, as e= 1 + h/2 + h^2/8, with numerical I(gtilde) and Q(gtilde)
    end
    
    for j=1:length(gam_range_z)
        Ez(1,j)= numericCOR(alp,bet,gam_range_z(j),gtil_choice(i)); 
        Ez(2,j)= 1 - gam_range_z(j)*C0 - gam_range_z(j)*gtil_choice(i)*C1 - C2*gam_range_z(j)^2; %O(gamma^2) COR with Schwager-Poschel
        esqr=1 - 2*bet*I*gam_range_z(j) + 2*(bet^2)*I*Q*(gam_range_z(j)^2);
        Ez(3,j)=sqrt(max([esqr,0])); %O(gamma^2) COR, as \sqrt(e^2), with numerical I(gtilde) and Q(gtilde)
        Ez(4,j)=1 - bet*I*gam_range_z(j) + (bet^2)*I*(Q - (1/2)*I)*(gam_range_z(j)^2); %O(gamma^2) COR, as e= 1 + h/2 + h^2/8, with numerical I(gtilde) and Q(gtilde)
        
        Az(1,j)= abs(Ez(2,j) - Ez(1,j)); %O(gamma^2) COR with Schwager-Poschel
        Az(2,j)= abs(Ez(3,j) - Ez(1,j)); %O(gamma^2) COR, as \sqrt(e^2), with numerical I(gtilde) and Q(gtilde)
        Az(3,j)= abs(Ez(4,j) - Ez(1,j)); %O(gamma^2) COR, as e= 1 + h/2 + h^2/8, with numerical I(gtilde) and Q(gtilde)
        
        Rz(1,j)= abs((Ez(2,j) - Ez(1,j))/Ez(1,j)); %O(gamma^2) COR with Schwager-Poschel
        Rz(2,j)= abs((Ez(3,j) - Ez(1,j))/Ez(1,j)); %O(gamma^2) COR, as \sqrt(e^2), with numerical I(gtilde) and Q(gtilde)
        Rz(3,j)= abs((Ez(4,j) - Ez(1,j))/Ez(1,j)); %O(gamma^2) COR, as e= 1 + h/2 + h^2/8, with numerical I(gtilde) and Q(gtilde)
    end
    
    fig=figure('Color','w');
    plot(gam_range,E(1,:),'-b','LineWidth',1);
    hold on
    plot(gam_range,E(2,:),'-r','LineWidth',1);
    hold on
    plot(gam_range,E(3,:),'--g','LineWidth',1);
    hold on
    plot(gam_range,E(4,:),'-.m','LineWidth',1);
    hold off
    title(strcat('$e$ w.r.t change in $\gamma$, with fixed $\tilde{g} = $',num2str(gtil_choice(i))),'Interpreter','latex','FontSize',20)
    legend('Numerical','SP - $O(\gamma^2)$, [Eq.~27]','$e=\sqrt{e^2}$ (Num. Coeffs.), [Eq.~91]','Linearized $e$ (Num. Coeffs.), [Eq.~92]','Interpreter','Latex','FontSize',16,'Location','Best')
    xlabel('$$\gamma$$', 'Interpreter','latex','FontSize',20)
    ylabel('$$e$$', 'Interpreter', 'latex', 'FontSize',20)
    
    pos=get(gca,'Position');
    scl=0.65; dims=[pos(3)-pos(1),pos(4)-pos(2)];
    ax_zoom=axes(fig,'OuterPosition',[pos(1),pos(2),scl*dims(1),scl*dims(2)],'Box','on');
    h(1)=plot(ax_zoom,gam_range_z,Ez(1,:),'-b','LineWidth',1);
    hold on
    h(2)=plot(ax_zoom,gam_range_z,Ez(2,:),'-r','LineWidth',1);
    hold on
    h(3)=plot(ax_zoom,gam_range_z,Ez(3,:),'--g','LineWidth',1);
    hold on 
    h(4)=plot(ax_zoom,gam_range_z,Ez(3,:),'-.m','LineWidth',1);
    hold off
    
    savefig(strcat(prefixA,'NumCompCOR',num2str(i),'.fig'));
    print(strcat(prefixB,'NumCompCOR',num2str(i),'.pdf'),'-dpdf');
    
%     fig=figure('Color','w');
%     plot(gam_range,R(1,:),'-r','LineWidth',1);
%     hold on 
%     plot(gam_range,R(2,:),'--g','LineWidth',1);
%     hold on 
%     plot(gam_range,R(3,:),'-.m','LineWidth',1);
%     hold off
%     title(strcat(' Relative errors for $\tilde{g} = $',num2str(gtil_choice(i))),'Interpreter','latex','FontSize',20)
%     legend('Schwager','$e=\sqrt{e^2}$ (Num. Coeffs.), [Eq.~91]','Linearized $e$ (Num. Coeffs.), [Eq.~92]','Interpreter','Latex','FontSize',16,'Location','Best')
%     xlabel('$$\gamma$$', 'Interpreter','latex','FontSize',20)
%     ylabel('Relative Error', 'Interpreter', 'latex', 'FontSize',20)
%     axis([gam_range(1),gam_range(end),0,min([10,max([R(1,:),R(2,:)])])]);
%     
%     pos=get(gca,'Position');
%     scl=0.65; dims=[pos(3)-pos(1),pos(4)-pos(2)];
%     ax_zoom=axes(fig,'OuterPosition',[pos(1),pos(2),scl*dims(1),scl*dims(2)],'Box','on');
%     h(1)=plot(ax_zoom,gam_range_z,Rz(1,:),'-r','LineWidth',1);
%     hold on
%     h(2)=plot(ax_zoom,gam_range_z,Rz(2,:),'--g','LineWidth',1);
%     hold on 
%     h(3)=plot(ax_zoom,gam_range_z,Rz(3,:),'-.m','LineWidth',1);
%     hold off
    
    fig=figure('Color','w');
    semilogy(gam_range,A(1,:),'-r','LineWidth',1);
    hold on 
    semilogy(gam_range,A(2,:),'--g','LineWidth',1);
    hold on 
    semilogy(gam_range,A(3,:),'-.m','LineWidth',1);
    hold off
    title(strcat(' Absolute errors for $\tilde{g} = $',num2str(gtil_choice(i))),'Interpreter','latex','FontSize',20)
    legend('SP - $O(\gamma^2)$, [Eq.~27]','$e=\sqrt{e^2}$ (Num. Coeffs.), [Eq.~91]','Linearized $e$ (Num. Coeffs.), [Eq.~92]','Interpreter','Latex','FontSize',16,'Location','Best')
    xlabel('$$\gamma$$', 'Interpreter','latex','FontSize',20)
    ylabel('Absolute Error', 'Interpreter', 'latex', 'FontSize',20)
    axis([gam_range(1),gam_range(end),0,min([10,max([A(1,:),A(2,:)])])]);
%     pos=get(gca,'Position');
%     scl=0.65; dims=[pos(3)-pos(1),pos(4)-pos(2)];
%     ax_zoom=axes(fig,'OuterPosition',[pos(1),pos(2),scl*dims(1),scl*dims(2)],'Box','on');
%     h(1)=plot(ax_zoom,gam_range_z,Az(1,:),'-r','LineWidth',1);
%     hold on
%     h(2)=plot(ax_zoom,gam_range_z,Az(2,:),'--g','LineWidth',1);
%     hold on 
%     h(3)=plot(ax_zoom,gam_range_z,Az(3,:),'-.m','LineWidth',1);
%     hold off
    
    savefig(strcat(prefixA,'NumCompCORErr',num2str(i),'.fig'));
    print(strcat(prefixB,'NumCompCORErr',num2str(i),'.pdf'),'-dpdf');
    
    fig=figure('Color','w');
    semilogy(gam_range,R(1,:),'-r','LineWidth',1);
    hold on 
    semilogy(gam_range,R(2,:),'--g','LineWidth',1);
    hold on 
    semilogy(gam_range,R(3,:),'-.m','LineWidth',1);
    hold off
    title(strcat(' Relative errors for $\tilde{g} = $',num2str(gtil_choice(i))),'Interpreter','latex','FontSize',20)
    legend('SP - $O(\gamma^2)$, [Eq.~27]','$e=\sqrt{e^2}$ (Num. Coeffs.), [Eq.~91]','Linearized $e$ (Num. Coeffs.), [Eq.~92]','Interpreter','Latex','FontSize',16,'Location','Best')
    xlabel('$$\gamma$$', 'Interpreter','latex','FontSize',20)
    ylabel('Relative Error', 'Interpreter', 'latex', 'FontSize',20)
    axis([gam_range(1),gam_range(end),0,min([10,max([R(1,:),R(2,:)])])]);
%     pos=get(gca,'Position');
%     scl=0.65; dims=[pos(3)-pos(1),pos(4)-pos(2)];
%     ax_zoom=axes(fig,'OuterPosition',[pos(1),pos(2),scl*dims(1),scl*dims(2)],'Box','on');
%     h(1)=plot(ax_zoom,gam_range_z,Rz(1,:),'-r','LineWidth',1);
%     hold on
%     h(2)=plot(ax_zoom,gam_range_z,Rz(2,:),'--g','LineWidth',1);
%     hold on 
%     h(3)=plot(ax_zoom,gam_range_z,Rz(3,:),'-.m','LineWidth',1);
%     hold off
    
    savefig(strcat(prefixA,'NumCompCORErrR',num2str(i),'.fig'));
    print(strcat(prefixB,'NumCompCORErrR',num2str(i),'.pdf'),'-dpdf');
end
end
