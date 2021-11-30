function CheckAnalVarCoeffCompare(choice) 
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


for i=1:length(gtil_choice)
    E=zeros(4,num_step);
    
    gam_range=linspace(0,min([gam_lim,1.1*max_gam(i)]),num_step);
    
    [In,Qn]=VarCORCoeffsNumeric(alp,bet,gtil_choice(i));
    [Ia,Qa]=VarCORCoeffsApprox(alp,bet,theta(i));
    
    if gtil_choice(i)==0
        In=((alp+1)/2)^(bet/(alp+1) - 1)*beta(3/2,bet/(alp+1));
        Qn=(1/2)*((alp+1)/2)^(bet/(alp+1) - 1)*beta(1/2,bet/(alp+1));
    end
    
    gam_cn=2/(bet*(Qn-(1/2)*In));
    gam_ca=2/(bet*(Qa-(1/2)*Ia));
    
    for j=1:length(gam_range)
        esqrn=1 - 2*bet*In*gam_range(j) + 2*(bet^2)*In*Qn*(gam_range(j)^2);
        elinn=1 - bet*In*gam_range(j) + (bet^2)*In*(Qn - (1/2)*In)*(gam_range(j)^2); 
        E(1,j) = min([sqrt(max([0,esqrn])), max([elinn,0])]);
        
        if gam_range(j)>gam_cn
            esqrn=1 - 2*bet*In*gam_range(j) + 2*(bet^2)*In*Qn*(gam_range(j)^2);
            E(2,j)=sqrt(max([0,esqrn])); %O(gamma^2) COR, as \sqrt(e^2), with numerical I(gtilde) and Q(gtilde)
        else
            E(2,j)=1 - bet*In*gam_range(j) + (bet^2)*In*(Qn - (1/2)*In)*(gam_range(j)^2);
        end
        
        esqra=1 - 2*bet*Ia*gam_range(j) + 2*(bet^2)*Ia*Qa*(gam_range(j)^2);
        elina=1 - bet*Ia*gam_range(j) + (bet^2)*Ia*(Qa - (1/2)*Ia)*(gam_range(j)^2);
        E(3,j) = min([ sqrt(max([0,esqra])), max([elina,0])]);
        
        if gam_range(j)>gam_ca
            esqra=1 - 2*bet*Ia*gam_range(j) + 2*(bet^2)*Ia*Qa*(gam_range(j)^2);
            E(4,j)=sqrt(max([0,esqra])); %O(gamma^2) COR, as \sqrt(e^2), with analytical I(gtilde) and Q(gtilde)
        else
            E(4,j)=1 - bet*Ia*gam_range(j) + (bet^2)*Ia*(Qa - (1/2)*Ia)*(gam_range(j)^2);
        end
    end
    
    figure('color','w')
    plot(gam_range,E(1,:)-E(2,:),'-b')
    hold on 
    plot(gam_range,E(3,:)-E(4,:),'--r')
    hold off 
    legend('Num. Coeffs.','Approx. Coeffs.','FontSize',12)
    title(strcat('Difference between the two threshold critera, $\alpha=',num2str(alp),'$, $\beta=',num2str(bet),...
        '$, $\tilde{g}=',num2str(gtil_choice(i)),'$'),'Interpreter','latex','FontSize',16);
    xlabel('$\gamma$','Interpreter','latex','FontSize',16)
    ylabel('$\Delta e = e_{new} - e_{old}$','Interpreter','latex','FontSize',16);
    
    

end

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
