function AnalVarCoeffCompare(choice) 
% Authors: Abhishek Chatterjee, Guillaume James, and Bernard Brogliato
% Address: Univ. Grenoble Alpes, INRIA, CNRS, Grenoble INP, LJK, Grenoble
%          38000 France 
%% Compare O(gamma^2) between Schwager-Poschel, numerical and analytical coeffeicients
alp=3/2;
if choice==1
    alp=3/2; bet=1; prefixA='Figures/matfig/AnalCoeff/case_1/'; prefixB='Figures/pdf/AnalCoeff/case_1/';
elseif choice==2
    alp=3/2; bet=5/4; prefixA='Figures/matfig/AnalCoeff/case_2/'; prefixB='Figures/pdf/AnalCoeff/case_2/';
elseif choice==3
    alp=3/2; bet=3/2;   prefixA='Figures/matfig/AnalCoeff/case_3/'; prefixB='Figures/pdf/AnalCoeff/case_3/';
elseif choice==4
    alp=3/2; bet=5/2;   prefixA='Figures/matfig/AnalCoeff/case_4/'; prefixB='Figures/pdf/AnalCoeff/case_4/';
elseif choice==5
    alp=1; bet=1; prefixA='Figures/matfig/AnalCoeff/case_5/'; prefixB='Figures/pdf/AnalCoeff/case_5/';
elseif choice==6
    alp=3/2; bet=5/4; prefixA='Figures/matfig/AnalCoeff/case_6/'; prefixB='Figures/pdf/AnalCoeff/case_6/';
else
    bet=3/2; prefixA='Figures/matfig/AnalCoeff/case_3/'; prefixB='Figures/pdf/AnalCoeff/case_3/';
end

% Get results for specific gtilde values 
gtil_choice=[0,0.05,1,10];
theta=[theta_gtilde(alp,gtil_choice(1)),theta_gtilde(alp,gtil_choice(2)),theta_gtilde(alp,gtil_choice(3)),theta_gtilde(alp,gtil_choice(4))];
%gtil_choice=[gtilde_theta(alp,theta(1)),gtilde_theta(alp,theta(2)),gtilde_theta(alp,theta(3)),gtilde_theta(alp,theta(4))];
%gtil_choice(1)=0;
z_fact=[3/1000,1/10,1/10000,1/1000];
max_gam=zeros(1,length(gtil_choice));
for i=1:length(gtil_choice)
    max_gam(i)=gamma_c(alp,bet,gtil_choice(i));
end

if choice==6
    gtil_choice=0;
end

num_step=1000;
gam_lim=0.1; 
if bet==3/2
    [C0,C1,C2]=ConstCORCoeffs(alp,bet);
else
    C2=0; % Order 2 constant coefficient is not known when bet~=3/2
    [C0,C1]=ConstCORCoeffs(alp,bet);
end

LowGamErr=[]; LowGamRng=[]; MidGamErr=[]; MidGamRng=[]; HighGamErr=[]; HighGamRng=[];
parts=[0.2,0.92];

tabletext={}; prec=3;
if bet==3/2
    ln=1; line="\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}"; tabletext{ln}=line;
    ln=ln+1; line="\\hline"; tabletext{ln} = line;
    ln=ln+1; line="& \\multicolumn{3}{c|}{Schwager-Poschel Coefficients} & \\multicolumn{3}{c|}{Numerical Coefficients} & \\multicolumn{3}{c|}{Analytical Coefficients} \\\\"; tabletext{ln}=line;
    ln=ln+1; line="\\hline"; tabletext{ln} = line;
    ln=ln+1; line=" $\\tilde{g}$ & Small $\\gamma$ & Medium $\\gamma$ & Large $\\gamma$  &  Small $\\gamma$ & Medium $\\gamma$ & Large $\\gamma$ & Small $\\gamma$ & Medium $\\gamma$ & Large $\\gamma$ \\\\"; tabletext{ln}=line;
    ln=ln+1; line="\\hline"; tabletext{ln} = line;
else
    ln=1; line="\\begin{tabular}{|c|c|c|c|c|c|c|}"; tabletext{ln}=line;
    ln=ln+1; line="\\hline"; tabletext{ln} = line;
    ln=ln+1; line="& \\multicolumn{3}{c|}{Numerical Coefficients} & \\multicolumn{3}{c|}{Analytical Coefficients} \\\\"; tabletext{ln}=line;
    ln=ln+1; line="\\hline"; tabletext{ln} = line;
    ln=ln+1; line=" $\\tilde{g}$ & Small $\\gamma$ & Medium $\\gamma$ & Large $\\gamma$  &  Small $\\gamma$ & Medium $\\gamma$ & Large $\\gamma$ \\\\"; tabletext{ln}=line;
    ln=ln+1; line="\\hline"; tabletext{ln} = line;
end


for i=1:length(gtil_choice)
    E=zeros(4,num_step); Ez=zeros(4,num_step);
    R=zeros(3,num_step); Rz=zeros(3,num_step);
    A=zeros(3,num_step); Az=zeros(3,num_step);
    
    gam_range=linspace(0,min([gam_lim,1.1*max_gam(i)]),num_step);
    gam_range_z=linspace(0,0.1*min([gam_lim,1.1*max_gam(i)]),num_step);
    
    [In,Qn]=VarCORCoeffsNumeric(alp,bet,gtil_choice(i));
    [Ia,Qa]=VarCORCoeffsApprox(alp,bet,theta(i));
    
    if gtil_choice(i)==0
        In=((alp+1)/2)^(bet/(alp+1) - 1)*beta(3/2,bet/(alp+1));
        Qn=(1/2)*((alp+1)/2)^(bet/(alp+1) - 1)*beta(1/2,bet/(alp+1));
    end
    
    gam_cn=2/(bet*(Qn-(1/2)*In));
    gam_ca=2/(bet*(Qa-(1/2)*Ia));
    
    for j=1:length(gam_range)
        E(1,j)= numericCOR(alp,bet,gam_range(j),gtil_choice(i));
        if choice==3
            E(2,j) = 1 - gam_range(j)*C0 - gam_range(j)*gtil_choice(i)*C1 - C2*gam_range(j)^2; %O(gamma^2) COR with Schwager-Poschel
        elseif choice==5
            E(2,j) = ExactCOR(1,gam_range(j),gtil_choice(i)); % Exact analytical solution of COR for alp=bet=1 
        elseif choice==6
            E(2,j) = ExactCOR(2,gam_range(j)); % Exact analytical solution of COR for alp=3/2 & bet=5/4, using Antypov-Elliott
        else
            E(2,j) = 0; % No analytical solution for other cases; 
        end
        
        esqrn=1 - 2*bet*In*gam_range(j) + 2*(bet^2)*In*Qn*(gam_range(j)^2);
        elinn=1 - bet*In*gam_range(j) + (bet^2)*In*(Qn - (1/2)*In)*(gam_range(j)^2); 
        E(3,j) = min([sqrt(max([0,esqrn])), max([elinn,0])]);
        
%         if gam_range(j)>gam_cn
%             esqrn=1 - 2*bet*In*gam_range(j) + 2*(bet^2)*In*Qn*(gam_range(j)^2);
%             E(3,j)=sqrt(max([0,esqrn])); %O(gamma^2) COR, as \sqrt(e^2), with numerical I(gtilde) and Q(gtilde)
%         else
%             E(3,j)=1 - bet*In*gam_range(j) + (bet^2)*In*(Qn - (1/2)*In)*(gam_range(j)^2);
%         end
        
        esqra=1 - 2*bet*Ia*gam_range(j) + 2*(bet^2)*Ia*Qa*(gam_range(j)^2);
        elina=1 - bet*Ia*gam_range(j) + (bet^2)*Ia*(Qa - (1/2)*Ia)*(gam_range(j)^2);
        E(4,j) = min([ sqrt(max([0,esqra])), max([elina,0])]);

%         if gam_range(j)>gam_ca
%             esqra=1 - 2*bet*Ia*gam_range(j) + 2*(bet^2)*Ia*Qa*(gam_range(j)^2);
%             E(4,j)=sqrt(max([0,esqra])); %O(gamma^2) COR, as \sqrt(e^2), with analytical I(gtilde) and Q(gtilde)
%         else
%             E(4,j)=1 - bet*Ia*gam_range(j) + (bet^2)*Ia*(Qa - (1/2)*Ia)*(gam_range(j)^2);
%         end
        
        
        if choice==5 || choice==6
            A(1,j)= abs(E(1,j) - E(2,j)); %O(gamma^2) COR of direct numerical, w.r.t. exact solution for the cases: 1. alp=bet=1, and 2. alp=3/2 & bet=5/4 
            A(2,j)= abs(E(3,j) - E(2,j)); %O(gamma^2) COR, as \sqrt(e^2), with numerical I(gtilde) and Q(gtilde), w.r.t. exact solution for the cases: 1. alp=bet=1, and 2. alp=3/2 & bet=5/4 
            A(3,j)= abs(E(4,j) - E(2,j)); %O(gamma^2) COR, as \sqrt(e^2), with analytical I(gtilde) and Q(gtilde), w.r.t. exact solution for the cases: 1. alp=bet=1, and 2. alp=3/2 & bet=5/4 
            
            R(1,j)= abs((E(1,j) - E(2,j))/E(2,j)); %O(gamma^2) COR using direct numerical integration, w.r.t. exact solution for the cases: 1. alp=bet=1, and 2. alp=3/2 & bet=5/4 
            R(2,j)= abs((E(3,j) - E(2,j))/E(2,j)); %O(gamma^2) COR, as \sqrt(e^2), with numerical I(gtilde) and Q(gtilde), w.r.t. exact solution for the cases: 1. alp=bet=1, and 2. alp=3/2 & bet=5/4 
            R(3,j)= abs((E(4,j) - E(2,j))/E(2,j)); %O(gamma^2) COR, as \sqrt(e^2), with analytical I(gtilde) and Q(gtilde), w.r.t. exact solution for the cases: 1. alp=bet=1, and 2. alp=3/2 & bet=5/4 
        else
            A(1,j)= abs(E(2,j) - E(1,j)); %O(gamma^2) COR with Schwager-Poschel
            A(2,j)= abs(E(3,j) - E(1,j)); %O(gamma^2) COR, as \sqrt(e^2), with numerical I(gtilde) and Q(gtilde)
            A(3,j)= abs(E(4,j) - E(1,j)); %O(gamma^2) COR, as \sqrt(e^2), with analytical I(gtilde) and Q(gtilde)
            
            R(1,j)= abs((E(2,j) - E(1,j))/E(1,j)); %O(gamma^2) COR with Schwager-Poschel
            R(2,j)= abs((E(3,j) - E(1,j))/E(1,j)); %O(gamma^2) COR, as \sqrt(e^2), with numerical I(gtilde) and Q(gtilde)
            R(3,j)= abs((E(4,j) - E(1,j))/E(1,j)); %O(gamma^2) COR, as \sqrt(e^2), with analytical I(gtilde) and Q(gtilde)
        end
    end
    
    for j=1:length(gam_range_z)
        Ez(1,j)= numericCOR(alp,bet,gam_range_z(j),gtil_choice(i));
        if choice==3
            Ez(2,j) = 1 - gam_range_z(j)*C0 - gam_range_z(j)*gtil_choice(i)*C1 - C2*gam_range_z(j)^2; %O(gamma^2) COR with Schwager-Poschel
        elseif choice==5
            Ez(2,j) = ExactCOR(1,gam_range_z(j),gtil_choice(i)); % Exact analytical solution of COR for alp=bet=1 
        elseif choice==6
            Ez(2,j) = ExactCOR(2,gam_range_z(j)); % Exact analytical solution of COR for alp=3/2 & bet=5/4, using Antypov-Elliott
        else
            Ez(2,j) = 0; % No analytical solution for other cases; 
        end
        
        esqrn=1 - 2*bet*In*gam_range_z(j) + 2*(bet^2)*In*Qn*(gam_range_z(j)^2);
        elinn=1 - bet*In*gam_range_z(j) + (bet^2)*In*(Qn - (1/2)*In)*(gam_range_z(j)^2);
        Ez(3,j) = min([sqrt(max([0,esqrn])), max([elinn,0])]);
        
%         if gam_range_z(j)>gam_cn
%             esqrn=1 - 2*bet*In*gam_range_z(j) + 2*(bet^2)*In*Qn*(gam_range_z(j)^2);
%             Ez(3,j)=sqrt(max([0,esqrn])); %O(gamma^2) COR, as \sqrt(e^2), with numerical I(gtilde) and Q(gtilde)
%         else
%             Ez(3,j)=1 - bet*In*gam_range_z(j) + (bet^2)*In*(Qn - (1/2)*In)*(gam_range_z(j)^2);
%         end
        esqra=1 - 2*bet*Ia*gam_range_z(j) + 2*(bet^2)*Ia*Qa*(gam_range_z(j)^2);
        elina=1 - bet*Ia*gam_range_z(j) + (bet^2)*Ia*(Qa - (1/2)*Ia)*(gam_range_z(j)^2);
        Ez(4,j) = min([ sqrt(max([0,esqra])), max([elina,0])]);
        
%         if gam_range_z(j)>gam_ca
%             esqra=1 - 2*bet*Ia*gam_range_z(j) + 2*(bet^2)*Ia*Qa*(gam_range_z(j)^2);
%             Ez(4,j)=sqrt(max([0,esqra])); %O(gamma^2) COR, as \sqrt(e^2), with analytical I(gtilde) and Q(gtilde)
%         else
%             Ez(4,j)=1 - bet*Ia*gam_range_z(j) + (bet^2)*Ia*(Qa - (1/2)*Ia)*(gam_range_z(j)^2);
%         end

        
        if choice==5 || choice==6
            Az(1,j)= abs(Ez(1,j) - Ez(2,j)); %O(gamma^2) COR of direct numerical, w.r.t. exact solution for the cases: 1. alp=bet=1, and 2. alp=3/2 & bet=5/4 
            Az(2,j)= abs(Ez(3,j) - Ez(2,j)); %O(gamma^2) COR, as \sqrt(e^2), with numerical I(gtilde) and Q(gtilde), w.r.t. exact solution for the cases: 1. alp=bet=1, and 2. alp=3/2 & bet=5/4 
            Az(3,j)= abs(Ez(4,j) - Ez(2,j)); %O(gamma^2) COR, as \sqrt(e^2), with analytical I(gtilde) and Q(gtilde), w.r.t. exact solution for the cases: 1. alp=bet=1, and 2. alp=3/2 & bet=5/4 
            
            Rz(1,j)= abs((Ez(1,j) - Ez(2,j))/Ez(2,j)); %O(gamma^2) COR using direct numerical integration, w.r.t. exact solution for the cases: 1. alp=bet=1, and 2. alp=3/2 & bet=5/4 
            Rz(2,j)= abs((Ez(3,j) - Ez(2,j))/Ez(2,j)); %O(gamma^2) COR, as \sqrt(e^2), with numerical I(gtilde) and Q(gtilde), w.r.t. exact solution for the cases: 1. alp=bet=1, and 2. alp=3/2 & bet=5/4 
            Rz(3,j)= abs((Ez(4,j) - Ez(2,j))/Ez(2,j)); %O(gamma^2) COR, as \sqrt(e^2), with analytical I(gtilde) and Q(gtilde), w.r.t. exact solution for the cases: 1. alp=bet=1, and 2. alp=3/2 & bet=5/4 
        else
            Az(1,j)= abs(Ez(2,j) - Ez(1,j)); %O(gamma^2) COR with Schwager-Poschel
            Az(2,j)= abs(Ez(3,j) - Ez(1,j)); %O(gamma^2) COR, as \sqrt(e^2), with numerical I(gtilde) and Q(gtilde)
            Az(3,j)= abs(Ez(4,j) - Ez(1,j)); %O(gamma^2) COR, as \sqrt(e^2), with analytical I(gtilde) and Q(gtilde)
            
            Rz(1,j)= abs((Ez(2,j) - Ez(1,j))/Ez(1,j)); %O(gamma^2) COR with Schwager-Poschel
            Rz(2,j)= abs((Ez(3,j) - Ez(1,j))/Ez(1,j)); %O(gamma^2) COR, as \sqrt(e^2), with numerical I(gtilde) and Q(gtilde)
            Rz(3,j)= abs((Ez(4,j) - Ez(1,j))/Ez(1,j)); %O(gamma^2) COR, as \sqrt(e^2), with analytical I(gtilde) and Q(gtilde)
        end
    end
    
    
    
    LowGamRng=[LowGamRng;[gam_range(1),parts(1)*gam_range(end)]]; MidGamRng=[MidGamRng;[parts(1)*gam_range(end),parts(2)*gam_range(end)]]; HighGamRng=[HighGamRng;[parts(2)*gam_range(end),gam_range(end)]];
    LowGamErr=[LowGamErr;[max(A(1,find(gam_range<parts(1)*gam_range(end)))),max(A(2,find(gam_range<parts(1)*gam_range(end)))), max(A(3,find(gam_range<parts(1)*gam_range(end)))) ]];
    MidGamErr=[MidGamErr; [ max(A(1,intersect(find(gam_range>=parts(1)*gam_range(end) ), find(gam_range<=parts(2)*gam_range(end))))), ...
        max(A(2,intersect(find(gam_range>=parts(1)*gam_range(end) ), find(gam_range<=parts(2)*gam_range(end))))), max(A(3,intersect(find(gam_range>=parts(1)*gam_range(end) ), find(gam_range<=parts(2)*gam_range(end)))))]];
    HighGamErr=[HighGamErr; [max(A(1,find(gam_range>parts(2)*gam_range(end)))), max(A(2,find(gam_range>parts(2)*gam_range(end)))), max(A(3,find(gam_range>parts(2)*gam_range(end))))]];
    
    
    if choice==3
        fig=figure('Color','w');
        plot(gam_range,E(1,:),'-b','LineWidth',1);
        hold on
        plot(gam_range,E(2,:),'--r','LineWidth',1);
        hold on
        plot(gam_range,E(3,:),'--g','LineWidth',1);
        hold on
        plot(gam_range,E(4,:),'-.m','LineWidth',1);
        hold off
        title(strcat('$e$ w.r.t change in $\gamma$, with fixed $\tilde{g} = ',num2str(gtil_choice(i)),'$, ($\alpha =',num2str(alp),'$ \& $\beta=',num2str(bet),'$)'),'Interpreter','latex','FontSize',16)
        if gtil_choice(i)==0
            legend('Num. Integ.,[Eq.~(4)]','Approx. (SP Coeff.) - ${\cal O}(\gamma^2)$, [Eq.~(25),(30)]','Using Eq.~(36) w/ Exact Coeffs., [Eq.~(39),(40)]','Using Eq.~(36) w/ Approx. Coeffs. [Eq.~(48),(49)]','Interpreter','Latex','FontSize',12,'Location','Best')
        else
            legend('Num. Integ.,[Eq.~(4)]','Approx. (SP Coeff.) - ${\cal O}(\gamma^2)$, [Eq.~(25),(30)]','Using Eq.~(36) w/ Num. Coeffs., [Eq.~(41),(42)]','Using Eq.~(36) w/ Approx. Coeffs. [Eq.~(48),(49)]','Interpreter','Latex','FontSize',12,'Location','Best')
        end
        xlabel('$$\gamma$$', 'Interpreter','latex','FontSize',20)
        ylabel('$$e$$', 'Interpreter', 'latex', 'FontSize',20)
        
        pos=get(gca,'Position');
        scl=0.65; dims=[pos(3)-pos(1),pos(4)-pos(2)];
        ax_zoom=axes(fig,'OuterPosition',[pos(1),pos(2),scl*dims(1),scl*dims(2)],'Box','on');
        h(1)=plot(ax_zoom,gam_range_z,Ez(1,:),'-b','LineWidth',1);
        hold on
        h(2)=plot(ax_zoom,gam_range_z,Ez(2,:),'--r','LineWidth',1);
        hold on
        h(3)=plot(ax_zoom,gam_range_z,Ez(3,:),'--g','LineWidth',1);
        hold on
        h(4)=plot(ax_zoom,gam_range_z,Ez(4,:),'-.m','LineWidth',1);
        hold off
        
        savefig(strcat(prefixA,'AnCompCOR',num2str(i),'.fig'));
        print(strcat(prefixB,'AnCompCOR',num2str(i),'.pdf'),'-dpdf');
        
        
        %---------------------------------------------------------%
        
        fig=figure('Color','w');
        semilogy(gam_range,A(1,:),'--r','LineWidth',1);
        hold on
        semilogy(gam_range,A(2,:),'--g','LineWidth',1);
        hold on
        semilogy(gam_range,A(3,:),'-.m','LineWidth',1);
        hold off
        title(strcat(' Absolute errors for $\tilde{g} = ',num2str(gtil_choice(i)),'$, ($\alpha =',num2str(alp),'$ \& $\beta=',num2str(bet),'$)'),'Interpreter','latex','FontSize',16)
        if gtil_choice(i)==0
            legend('Approx. (SP Coeff.) - ${\cal O}(\gamma^2)$, [Eq.~(25),(30)]','Using Eq.~(36) w/ Exact Coeffs., [Eq.~(39),(40)]','Using Eq.~(36) w/ Approx. Coeffs. [Eq.~(48),(49)]','Interpreter','Latex','FontSize',12,'Location','Best')
        else
            legend('Approx. (SP Coeff.) - ${\cal O}(\gamma^2)$, [Eq.~(25),(30)]','Using Eq.~(36) w/ Num. Coeffs., [Eq.~(41),(42)]','Using Eq.~(36) w/ Approx. Coeffs. [Eq.~(48),(49)]','Interpreter','Latex','FontSize',12,'Location','Best')
        end
        xlabel('$$\gamma$$', 'Interpreter','latex','FontSize',20)
        ylabel('Absolute Error', 'Interpreter', 'latex', 'FontSize',20)
        axis([gam_range(1),gam_range(end),0,min([100,max([A(1,:),A(2,:)])])]);
        
%         pos=get(gca,'Position');
%         scl=0.65; dims=[pos(3)-pos(1),pos(4)-pos(2)];
%         ax_zoom=axes(fig,'OuterPosition',[pos(1),pos(2),scl*dims(1),scl*dims(2)],'Box','on');
%         h(1)=plot(ax_zoom,gam_range_z,Az(1,:),'--r','LineWidth',1);
%         hold on
%         h(2)=plot(ax_zoom,gam_range_z,Az(2,:),'--g','LineWidth',1);
%         hold on
%         h(3)=plot(ax_zoom,gam_range_z,Az(3,:),'-.m','LineWidth',1);
%         hold off
        
        
        savefig(strcat(prefixA,'AnCompCORErr',num2str(i),'.fig'));
        print(strcat(prefixB,'AnCompCORErr',num2str(i),'.pdf'),'-dpdf');
        
        %----------------------------------------------------------%
        
        fig=figure('Color','w');
        semilogy(gam_range,R(1,:),'--r','LineWidth',1);
        hold on
        semilogy(gam_range,R(2,:),'--g','LineWidth',1);
        hold on
        semilogy(gam_range,R(3,:),'-.m','LineWidth',1);
        hold off
        title(strcat(' Relative errors for $\tilde{g} = ',num2str(gtil_choice(i)),'$, ($\alpha =',num2str(alp),'$ \& $\beta=',num2str(bet),'$)'),'Interpreter','latex','FontSize',16)
        if gtil_choice(i)==0
            legend('Approx. (SP Coeff.) - ${\cal O}(\gamma^2)$, [Eq.~(25),(30)]','Using Eq.~(36) w/ Exact Coeffs., [Eq.~(39),(40)]','Using Eq.~(36) w/ Approx. Coeffs. [Eq.~(48),(49)]','Interpreter','Latex','FontSize',12,'Location','Best')
        else
            legend('Approx. (SP Coeff.) - ${\cal O}(\gamma^2)$, [Eq.~(25),(30)]','Using Eq.~(36) w/ Num. Coeffs., [Eq.~(41),(42)]','Using Eq.~(36) w/ Approx. Coeffs. [Eq.~(48),(49)]','Interpreter','Latex','FontSize',12,'Location','Best')
        end
        xlabel('$$\gamma$$', 'Interpreter','latex','FontSize',20)
        ylabel('Relative Error', 'Interpreter', 'latex', 'FontSize',20)
        axis([gam_range(1),gam_range(end),0,min([100,max([R(1,:),R(2,:)])])]);
        
%         pos=get(gca,'Position');
%         scl=0.65; dims=[pos(3)-pos(1),pos(4)-pos(2)];
%         ax_zoom=axes(fig,'OuterPosition',[pos(1),pos(2),scl*dims(1),scl*dims(2)],'Box','on');
%         h(1)=plot(ax_zoom,gam_range_z,Rz(1,:),'--r','LineWidth',1);
%         hold on
%         h(2)=plot(ax_zoom,gam_range_z,Rz(2,:),'--g','LineWidth',1);
%         hold on
%         h(3)=plot(ax_zoom,gam_range_z,Rz(3,:),'-.m','LineWidth',1);
%         hold off
        
        savefig(strcat(prefixA,'AnCompCORErrR',num2str(i),'.fig'));
        print(strcat(prefixB,'AnCompCORErrR',num2str(i),'.pdf'),'-dpdf');
        
        %-----------------------------------------------------------%
        
        ln=ln+1;
        line = strcat( "\\multirow{2}{*}{$", textstr(gtil_choice(i),prec), "$}"," & $\\gamma < ", textstr(LowGamRng(i,2),prec), "$ & $", textstr(MidGamRng(i,1),prec) , " \\leq \\gamma \\leq ", textstr(MidGamRng(i,2),prec) ...
            , "$ & $ \\gamma > " , textstr(HighGamRng(i,1),prec), " $ & $ \\gamma < ", textstr(LowGamRng(i,2),prec), " $ & $ ", textstr(MidGamRng(i,1),prec) , " \\leq \\gamma \\leq ", textstr(MidGamRng(i,2),prec), "$ & $ \\gamma > " , textstr(HighGamRng(1),prec)  ...
            ," $ & $ \\gamma < ", textstr(LowGamRng(i,2),prec), " $ & $ ", textstr(MidGamRng(i,1),prec) , " \\leq \\gamma \\leq ", textstr(MidGamRng(i,2),prec), " $ & $ \\gamma > " , textstr(HighGamRng(i,1),prec), "$ \\\\" );
        tabletext{ln} = line;
        ln=ln+1;
        line="\\cline{2-10}";
        tabletext{ln}=line;
        ln=ln+1;
        line = strcat( " & $ ", textstr(LowGamErr(1),prec), " $ & $ ", textstr(MidGamErr(1),prec), " $ & $ ", textstr(HighGamErr(1),prec), "$ & $",  textstr(LowGamErr(2),prec), " $ & $ ", textstr(MidGamErr(2),prec), " $ & $ ", textstr(HighGamErr(2),prec),...
            "$ & $", textstr(LowGamErr(3),prec), "$ & $", textstr(MidGamErr(3),prec), "$ & $", textstr(HighGamErr(3),prec), '$ \\\\');
        tabletext{ln} = line;
        ln=ln+1;
        line="\\hline";
        tabletext{ln} = line;
        
    elseif choice==5 || choice==6
        fig=figure('Color','w');
        plot(gam_range,E(2,:),'-b','LineWidth',1);
        hold on
        plot(gam_range,E(1,:),'--r','LineWidth',1);
        hold on
        plot(gam_range,E(3,:),'--g','LineWidth',1);
        hold on
        plot(gam_range,E(4,:),'-.m','LineWidth',1);
        hold off
        title(strcat('$e$ w.r.t change in $\gamma$, with fixed $\tilde{g} = ',num2str(gtil_choice(i)),'$, ($\alpha =',num2str(alp),'$ \& $\beta=',num2str(bet),'$)'),'Interpreter','latex','FontSize',16)
        if gtil_choice(i)==0
            legend('Known Analytic. Sol.','Num. Integ.,[Eq.~(4)]','Using Eq.~(36) w/ Exact Coeffs., [Eq.~(39),(40)]','Using Eq.~(36) w/ Approx. Coeffs. [Eq.~(48),(49)]','Interpreter','Latex','FontSize',12,'Location','Best')
        else
            legend('Known Analytic. Sol.','Num. Integ.,[Eq.~(4)]','Using Eq.~(36) w/ Num. Coeffs., [Eq.~(41),(42)]','Using Eq.~(36) w/ Approx. Coeffs. [Eq.~(48),(49)]','Interpreter','Latex','FontSize',12,'Location','Best')
        end
        xlabel('$$\gamma$$', 'Interpreter','latex','FontSize',20)
        ylabel('$$e$$', 'Interpreter', 'latex', 'FontSize',20)
        
%         pos=get(gca,'Position');
%         scl=0.65; dims=[pos(3)-pos(1),pos(4)-pos(2)];
%         ax_zoom=axes(fig,'OuterPosition',[pos(1),pos(2),scl*dims(1),scl*dims(2)],'Box','on');
%         h(1)=plot(ax_zoom,gam_range_z,Ez(2,:),'-b','LineWidth',1);
%         hold on
%         h(2)=plot(ax_zoom,gam_range_z,Ez(1,:),'--r','LineWidth',1);
%         hold on
%         h(3)=plot(ax_zoom,gam_range_z,Ez(3,:),'--g','LineWidth',1);
%         hold on
%         h(4)=plot(ax_zoom,gam_range_z,Ez(4,:),'-.m','LineWidth',1);
%         hold off
        
        savefig(strcat(prefixA,'AnCompCOR',num2str(i),'.fig'));
        print(strcat(prefixB,'AnCompCOR',num2str(i),'.pdf'),'-dpdf');
        
        %--------------------------------------------------------------%
        
        fig=figure('Color','w');
        semilogy(gam_range,A(1,:),'--r','LineWidth',1);
        hold on
        semilogy(gam_range,A(2,:),'--g','LineWidth',1);
        hold on
        semilogy(gam_range,A(3,:),'-.m','LineWidth',1);
        hold off
        title(strcat(' Absolute errors for $\tilde{g} = ',num2str(gtil_choice(i)),'$, ($\alpha =',num2str(alp),'$ \& $\beta=',num2str(bet),'$)'),'Interpreter','latex','FontSize',16)
        if gtil_choice(i)==0
            legend('Num. Integ.,[Eq.~(4)]','Using Eq.~(36) w/ Exact Coeffs., [Eq.~(39),(40)]','Using Eq.~(36) w/ Approx. Coeffs. [Eq.~(48),(49)]','Interpreter','Latex','FontSize',12,'Location','Best')
        else
            legend('Num. Integ.,[Eq.~(4)]','Using Eq.~(36) w/ Num. Coeffs., [Eq.~(41),(42)]','Using Eq.~(36) w/ Approx. Coeffs. [Eq.~(48),(49)]','Interpreter','Latex','FontSize',12,'Location','Best')
        end
        xlabel('$$\gamma$$', 'Interpreter','latex','FontSize',20)
        ylabel('Absolute Error', 'Interpreter', 'latex', 'FontSize',20)
        axis([gam_range(1),gam_range(end),0,min([100,max([A(1,:),A(2,:)])])]);
        
%         pos=get(gca,'Position');
%         scl=0.65; dims=[pos(3)-pos(1),pos(4)-pos(2)];
%         ax_zoom=axes(fig,'OuterPosition',[pos(1),pos(2),scl*dims(1),scl*dims(2)],'Box','on');
%         h(1)=plot(ax_zoom,gam_range_z,Az(1,:),'--r','LineWidth',1);
%         hold on
%         h(2)=plot(ax_zoom,gam_range_z,Az(2,:),'--g','LineWidth',1);
%         hold on
%         h(3)=plot(ax_zoom,gam_range_z,Az(3,:),'-.m','LineWidth',1);
%         hold off
        
        savefig(strcat(prefixA,'AnCompCORErr',num2str(i),'.fig'));
        print(strcat(prefixB,'AnCompCORErr',num2str(i),'.pdf'),'-dpdf');figure('Color','w');
        
        %--------------------------------------------------------------%
        
        fig=figure('Color','w');
        semilogy(gam_range,R(1,:),'--r','LineWidth',1);
        hold on
        semilogy(gam_range,R(2,:),'--g','LineWidth',1);
        hold on
        semilogy(gam_range,R(3,:),'-.m','LineWidth',1);
        hold off
        title(strcat(' Relative errors for $\tilde{g} = ',num2str(gtil_choice(i)),'$, ($\alpha =',num2str(alp),'$ \& $\beta=',num2str(bet),'$)'),'Interpreter','latex','FontSize',16)
        if gtil_choice(i)==0
            legend('Num. Integ.,[Eq.~(4)]','Using Eq.~(36) w/ Exact Coeffs., [Eq.~(39),(40)]','Using Eq.~(36) w/ Approx. Coeffs. [Eq.~(48),(49)]','Interpreter','Latex','FontSize',12,'Location','Best')
        else
            legend('Num. Integ.,[Eq.~(4)]','Using Eq.~(36) w/ Num. Coeffs., [Eq.~(41),(42)]','Using Eq.~(36) w/ Approx. Coeffs. [Eq.~(48),(49)]','Interpreter','Latex','FontSize',12,'Location','Best')
        end
        xlabel('$$\gamma$$', 'Interpreter','latex','FontSize',20)
        ylabel('Relative Error', 'Interpreter', 'latex', 'FontSize',20)
        axis([gam_range(1),gam_range(end),0,min([100,max([R(1,:),R(2,:)])])]);
        
%         pos=get(gca,'Position');
%         scl=0.65; dims=[pos(3)-pos(1),pos(4)-pos(2)];
%         ax_zoom=axes(fig,'OuterPosition',[pos(1),pos(2),scl*dims(1),scl*dims(2)],'Box','on');
%         h(1)=plot(ax_zoom,gam_range_z,Rz(1,:),'--r','LineWidth',1);
%         hold on
%         h(2)=plot(ax_zoom,gam_range_z,Rz(2,:),'--g','LineWidth',1);
%         hold on
%         h(3)=plot(ax_zoom,gam_range_z,Rz(3,:),'-.m','LineWidth',1);
%         hold off
        
        savefig(strcat(prefixA,'AnCompCORErrR',num2str(i),'.fig'));
        print(strcat(prefixB,'AnCompCORErrR',num2str(i),'.pdf'),'-dpdf');
        
        %--------------------------------------------------------------%

        ln=ln+1;
        line = strcat( "\\multirow{2}{*}{$", textstr(gtil_choice(i),prec), "$}"," & $\\gamma < ", textstr(LowGamRng(i,2),prec), "$ & $", textstr(MidGamRng(i,1),prec) , " \\leq \\gamma \\leq ", textstr(MidGamRng(i,2),prec) ...
            , "$ & $ \\gamma > " , textstr(HighGamRng(i,1),prec), " $ & $ \\gamma < ", textstr(LowGamRng(i,2),prec), " $ & $ ", textstr(MidGamRng(i,1),prec) , " \\leq \\gamma \\leq ", textstr(MidGamRng(i,2),prec), "$ & $ \\gamma > " , textstr(HighGamRng(1),prec)  ...
            ," $ & $ \\gamma < ", textstr(LowGamRng(i,2),prec), " $ & $ ", textstr(MidGamRng(i,1),prec) , " \\leq \\gamma \\leq ", textstr(MidGamRng(i,2),prec), " $ & $ \\gamma > " , textstr(HighGamRng(i,1),prec), "$ \\\\" );
        tabletext{ln} = line;
        ln=ln+1;
        line="\\cline{2-10}";
        tabletext{ln}=line;
        ln=ln+1;
        line = strcat( " & $ ", textstr(LowGamErr(1),prec), " $ & $ ", textstr(MidGamErr(1),prec), " $ & $ ", textstr(HighGamErr(1),prec), "$ & $",  textstr(LowGamErr(2),prec), " $ & $ ", textstr(MidGamErr(2),prec), " $ & $ ", textstr(HighGamErr(2),prec),...
            "$ & $", textstr(LowGamErr(3),prec), "$ & $", textstr(MidGamErr(3),prec), "$ & $", textstr(HighGamErr(3),prec), '$ \\\\');
        tabletext{ln} = line;
        ln=ln+1;
        line="\\hline";
        tabletext{ln} = line;
    else
        fig=figure('Color','w');
        plot(gam_range,E(1,:),'-b','LineWidth',1);
        hold on
        plot(gam_range,E(3,:),'--g','LineWidth',1);
        hold on
        plot(gam_range,E(4,:),'-.m','LineWidth',1);
        hold off
        title(strcat('$e$ w.r.t change in $\gamma$, with fixed $\tilde{g} =',num2str(gtil_choice(i)),'$, ($\alpha =',num2str(alp),'$ \& $\beta=',num2str(bet),'$)'),'Interpreter','latex','FontSize',16)
        if gtil_choice(i)==0
            legend('Num. Integ.,[Eq.~(4)]','Using Eq.~(36) w/ Exact Coeffs., [Eq.~(39),(40)]','Using Eq.~(36) w/ Approx. Coeffs. [Eq.~(48),(49)]','Interpreter','Latex','FontSize',12,'Location','Best')
        else
            legend('Num. Integ.,[Eq.~(4)]','Using Eq.~(36) w/ Num. Coeffs., [Eq.~(41),(42)]','Using Eq.~(36) w/ Approx. Coeffs. [Eq.~(48),(49)]','Interpreter','Latex','FontSize',12,'Location','Best')
        end
        xlabel('$$\gamma$$', 'Interpreter','latex','FontSize',20)
        ylabel('$$e$$', 'Interpreter', 'latex', 'FontSize',20)
        
        pos=get(gca,'Position');
        scl=0.65; dims=[pos(3)-pos(1),pos(4)-pos(2)];
        ax_zoom=axes(fig,'OuterPosition',[pos(1),pos(2),scl*dims(1),scl*dims(2)],'Box','on');
        h(1)=plot(ax_zoom,gam_range_z,Ez(1,:),'-b','LineWidth',1);
        hold on
        h(2)=plot(ax_zoom,gam_range_z,Ez(3,:),'--g','LineWidth',1);
        hold on
        h(3)=plot(ax_zoom,gam_range_z,Ez(4,:),'-.m','LineWidth',1);
        hold off
        
        savefig(strcat(prefixA,'AnCompCOR',num2str(i),'.fig'));
        print(strcat(prefixB,'AnCompCOR',num2str(i),'.pdf'),'-dpdf');
        
        %--------------------------------------------------------------%
        
        fig=figure('Color','w');
        semilogy(gam_range,A(2,:),'--g','LineWidth',1);
        hold on
        semilogy(gam_range,A(3,:),'-.m','LineWidth',1);
        hold off
        title(strcat(' Absolute errors for $\tilde{g} = ',num2str(gtil_choice(i)),'$, ($\alpha =',num2str(alp),'$ \& $\beta=',num2str(bet),'$)'),'Interpreter','latex','FontSize',16)
        if gtil_choice(i)==0
            legend('Using Eq.~(36) w/ Exact Coeffs., [Eq.~(39),(40)]','Using Eq.~(36) w/ Approx. Coeffs. [Eq.~(48),(49)]','Interpreter','Latex','FontSize',12,'Location','Best')
        else
            legend('Using Eq.~(36) w/ Num. Coeffs., [Eq.~(41),(42)]','Using Eq.~(36) w/ Approx. Coeffs. [Eq.~(48),(49)]','Interpreter','Latex','FontSize',12,'Location','Best')
        end
        xlabel('$$\gamma$$', 'Interpreter','latex','FontSize',20)
        ylabel('Absolute Error', 'Interpreter', 'latex', 'FontSize',20)
        axis([gam_range(1),gam_range(end),0,min([100,max([A(1,:),A(2,:)])])]);
        
%         pos=get(gca,'Position');
%         scl=0.65; dims=[pos(3)-pos(1),pos(4)-pos(2)];
%         ax_zoom=axes(fig,'OuterPosition',[pos(1),pos(2),scl*dims(1),scl*dims(2)],'Box','on');
%         h(1)=plot(ax_zoom,gam_range_z,Az(2,:),'--g','LineWidth',1);
%         hold on
%         h(2)=plot(ax_zoom,gam_range_z,Az(3,:),'-.m','LineWidth',1);
%         hold off
        
        savefig(strcat(prefixA,'AnCompCORErr',num2str(i),'.fig'));
        print(strcat(prefixB,'AnCompCORErr',num2str(i),'.pdf'),'-dpdf');
        
        %--------------------------------------------------------------%
        
        
        fig=figure('Color','w');
        semilogy(gam_range,R(2,:),'--g','LineWidth',1);
        hold on
        semilogy(gam_range,R(3,:),'-.m','LineWidth',1);
        hold off
        title(strcat(' Relative errors for $\tilde{g} = ',num2str(gtil_choice(i)),'$, ($\alpha =',num2str(alp),'$ \& $\beta=',num2str(bet),'$)'),'Interpreter','latex','FontSize',16)
        if gtil_choice(i)==0
            legend('Using Eq.~(36) w/ Exact Coeffs., [Eq.~(39),(40)]','Using Eq.~(36) w/ Approx. Coeffs. [Eq.~(48),(49)]','Interpreter','Latex','FontSize',12,'Location','Best')
        else
            legend('Using Eq.~(36) w/ Num. Coeffs., [Eq.~(41),(42)]','Using Eq.~(36) w/ Approx. Coeffs. [Eq.~(48),(49)]','Interpreter','Latex','FontSize',12,'Location','Best')
        end
        xlabel('$$\gamma$$', 'Interpreter','latex','FontSize',20)
        ylabel('Relative Error', 'Interpreter', 'latex', 'FontSize',20)
        axis([gam_range(1),gam_range(end),0,min([100,max([R(1,:),R(2,:)])])]);
        
%         pos=get(gca,'Position');
%         scl=0.65; dims=[pos(3)-pos(1),pos(4)-pos(2)];
%         ax_zoom=axes(fig,'OuterPosition',[pos(1),pos(2),scl*dims(1),scl*dims(2)],'Box','on');
%         h(1)=plot(ax_zoom,gam_range_z,Rz(2,:),'--g','LineWidth',1);
%         hold on
%         h(2)=plot(ax_zoom,gam_range_z,Rz(1,:),'-.m','LineWidth',1);
%         hold off
        
        savefig(strcat(prefixA,'AnCompCORErrR',num2str(i),'.fig'));
        print(strcat(prefixB,'AnCompCORErrR',num2str(i),'.pdf'),'-dpdf');
        
        %--------------------------------------------------------------%
        
        ln=ln+1;
        line = strcat( "\\multirow{2}{*}{$", textstr(gtil_choice(i),prec), "$}"," & $ \\gamma < ", textstr(LowGamRng(i,2),prec), " $ & $ ", textstr(MidGamRng(i,1),prec) , " \\leq \\gamma \\leq ", textstr(MidGamRng(i,2),prec), "$ & $ \\gamma > " , textstr(HighGamRng(i,1),prec)  ...
            ," $ & $ \\gamma < ", textstr(LowGamRng(i,2),prec), " $ & $ ", textstr(MidGamRng(i,1),prec) , " \\leq \\gamma \\leq ", textstr(MidGamRng(i,2),prec), " $ & $ \\gamma > " , textstr(HighGamRng(i,1),prec), "$ \\\\" );
        tabletext{ln} = line;
        ln=ln+1;
        line="\\cline{2-7}";
        tabletext{ln}=line;
        ln=ln+1;
        line = strcat( " & $ ", textstr(LowGamErr(2),prec), " $ & $ ", textstr(MidGamErr(2),prec), " $ & $ ", textstr(HighGamErr(2),prec),...
            "$ & $", textstr(LowGamErr(3),prec), "$ & $", textstr(MidGamErr(3),prec), "$ & $", textstr(HighGamErr(3),prec), '$ \\\\');
        tabletext{ln} = line;
        ln=ln+1;
        line="\\hline";
        tabletext{ln} = line;
    end
    
    

end

ln=ln+1; line="\\end{tabular}"; tabletext{ln} = line;
ln=ln+1; line=strcat("\\caption{ Maximum Absoulute Errors for low, medium and high $\\gamma$ values, given $\\alpha=",num2str(alp),"$, and $\\beta = ",num2str(bet),"$}"); tabletext{ln} = line;

for i=1:length(tabletext)
    fprintf(tabletext{i}); 
    fprintf("\n")
end
fprintf("\n")
fprintf("\n")
fprintf("\n")
end

function vstr=textstr(val,p)
vstr=num2str(val,p); 
if contains(vstr,'e')
    C=strsplit(vstr,'e');
    vstr=strcat(C{1},' \\times 10^{',C{2},'}');
end
end

function e=ExactCOR(choice,gam,varargin)
gtil=0;
if choice==1
    if ~isempty(varargin)
        gtil=varargin{1};
    end
    if gtil==0
        xi=gam/2; wn = sqrt(4-gam^2)/2;
        e=exp(-xi*pi/wn);
    else
        [tmx,tmn] = linstationarytime(gam,gtil);
        [umn,~]=linsol(gam,gtil,tmn);
        if umn<=0
            tc = 2*tmx;
            [uc,duc] = linsol(gam,gtil,tc);
            rho=0.1; maxiter=100; iter=1;
            while abs(uc)>1e-4 && abs(duc)>1e-4 && iter<=maxiter 
                tc = tc - rho*uc/duc;
                [uc,duc] = linsol(gam,gtil,tc);
                iter=iter+1;
            end
            e = - duc;
        else
            e=0;
        end
    end
elseif choice==2
    tc=mapcoltime(gam);
    [~,duf] = mapsol(gam,tc);
    e = - duf;
else
    error('Invalid Choice')
end
end
function [u,du] = linsol(gam,gtil,t)
xi=gam/2; wn = sqrt(4-gam^2)/2;
u = -gtil*exp(-xi*t)*cos(wn*t) + ((1-xi*gtil)/wn)*exp(-xi*t)*sin(wn*t) + gtil;
du = exp(-xi*t)*cos(wn*t) + ( wn*gtil - (xi*(1-xi*gtil))/(wn) )*exp(-xi*t)*sin(wn*t);
end
function [tmax,tmin] = linstationarytime(gam,gtil)
xi=gam/2; wn = sqrt(4-gam^2)/2;
% tmax = (1/wn)*(atan((xi*(1-xi*gtil) - wn^2*gtil)/wn) + pi);
% tmin = (1/wn)*(atan((xi*(1-xi*gtil) - wn^2*gtil)/wn) + 2*pi);
tmax = (1/wn)*(atan((wn^2*gtil - xi*(1-xi*gtil))/wn) + pi/2);
tmin = (1/wn)*(atan((wn^2*gtil - xi*(1-xi*gtil))/wn) + 3*pi/2);
end
function [u,du]=mapsol(gam,t)
A=(5/4)^(2/5); n=4/5;
[ul,dul]=linsol((sqrt(5)/2)*gam,0,t);
u=A*(ul^n);
du=dul;
%du=n*A*(ul^(n-1))*dul;
end
function tc=mapcoltime(gam)
wn = sqrt(4-(5/4)*gam^2)/2;
tc=pi/wn;
end
