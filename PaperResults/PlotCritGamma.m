function PlotCritGamma()
addpath ../
%{
alp=1; bet=1:0.01:1.5;
gam_c=zeros(2,length(bet));
for i =1:length(bet)
    gam_c(1,i)=SearchCritGamma(alp,bet(i),1e-4,1);
    gam_c(2,i)=findcritgamma_0g(alp,bet(i),1e-4,0.1);
    fprintf('beta = %.4f\n',bet(i))
end
bet_k=alp; gam_k=2*(bet_k^((-3/(2*(alp+1))) - 1));  %gam_k=sqrt(4/((3/2)^(2/5)*((2/3)^(4/5))));
bet_a=(alp+1)/2; gam_a=2*(bet_a^(1-bet_a/(alp+1)));
plot(bet,gam_c(1,:),'-b',bet,gam_c(2,:),'-r',bet_k,gam_k,'og',bet_a,gam_a,'om');
title(strcat('Critical $\gamma_c$ for $\alpha =',num2str(alp),'$'),'FontSize',20,'Interpreter','latex')
xlabel('$\beta$','FontSize',20,'Interpreter','latex')
ylabel('$\gamma_c$','FontSize',20,'Interpreter','latex')
legend('Binary','Secant','Known Unscaled','Known Scaled (Antypov)','FontSize',20,'Interpreter','latex')
print('TempFigures/critgamma1.pdf', '-dpdf');

alp=5/4; bet=1:0.01:1.5;
gam_c=zeros(2,length(bet));
for i =1:length(bet)
    gam_c(1,i)=SearchCritGamma(alp,bet(i),1e-4,1);
    gam_c(2,i)=findcritgamma_0g(alp,bet(i),1e-4,0.1);
    fprintf('beta = %.4f\n',bet(i))
end
bet_k=alp; gam_k=2*(bet_k^((-3/(2*(alp+1))) - 1));  %gam_k=sqrt(4/((3/2)^(2/5)*((2/3)^(4/5))));
bet_a=(alp+1)/2; gam_a=2*(bet_a^(1-bet_a/(alp+1)));
plot(bet,gam_c(1,:),'-b',bet,gam_c(2,:),'-r',bet_k,gam_k,'og',bet_a,gam_a,'om');
title(strcat('Critical $\gamma_c$ for $\alpha =',num2str(alp),'$'),'FontSize',20,'Interpreter','latex')
xlabel('$\beta$','FontSize',20,'Interpreter','latex')
ylabel('$\gamma_c$','FontSize',20,'Interpreter','latex')
legend('Binary','Secant','Known Unscaled','Known Scaled (Antypov)','FontSize',20,'Interpreter','latex')
print('TempFigures/critgamma2.pdf', '-dpdf');
%}
alp=3/2; bet=1:0.01:1.5;
gam_c=zeros(2,length(bet));
for i =1:length(bet)
    gam_c(1,i)=SearchCritGamma(alp,bet(i),1e-4,1);
    gam_c(2,i)=findcritgamma_0g(alp,bet(i),1e-4,0.1);
    fprintf('beta = %.4f\n',bet(i))
end
bet_k=alp; gam_k=2*(bet_k^((-3/(2*(alp+1))) - 1));  %gam_k=sqrt(4/((3/2)^(2/5)*((2/3)^(4/5))));
bet_a=(alp+1)/2; gam_a=2*(bet_a^(1-bet_a/(alp+1)));
plot(bet,gam_c(1,:),'-b',bet,gam_c(2,:),'-r',bet_k,gam_k,'og',bet_a,gam_a,'om');
title(strcat('Critical $\gamma_c$ for $\alpha =',num2str(alp),'$'),'FontSize',20,'Interpreter','latex')
xlabel('$\beta$','FontSize',20,'Interpreter','latex')
ylabel('$\gamma_c$','FontSize',20,'Interpreter','latex')
legend('Binary','Secant','Known Unscaled','Known Scaled (Antypov)','FontSize',20,'Interpreter','latex')
print('TempFigures/critgamma3.pdf', '-dpdf');
%{
alp=5/2; bet=1:0.01:1.5;
gam_c=zeros(2,length(bet));
for i =1:length(bet)
    gam_c(1,i)=SearchCritGamma(alp,bet(i),1e-4,1);
    gam_c(2,i)=findcritgamma_0g(alp,bet(i),1e-4,0.1);
    fprintf('beta = %.4f\n',bet(i))
end
bet_k=alp; gam_k=2*(bet_k^((-3/(2*(alp+1))) - 1));  %gam_k=sqrt(4/((3/2)^(2/5)*((2/3)^(4/5))));
bet_a=(alp+1)/2; gam_a=2*(bet_a^(1-bet_a/(alp+1)));
plot(bet,gam_c(1,:),'-b',bet,gam_c(2,:),'-r',bet_k,gam_k,'og',bet_a,gam_a,'om');
title(strcat('Critical $\gamma_c$ for $\alpha =',num2str(alp),'$'),'FontSize',20,'Interpreter','latex')
xlabel('$\beta$','FontSize',20,'Interpreter','latex')
ylabel('$\gamma_c$','FontSize',20,'Interpreter','latex')
legend('Binary','Secant','Known Unscaled','Known Scaled (Antypov)','FontSize',20,'Interpreter','latex')
print('TempFigures/critgamma4.pdf', '-dpdf');

%}
rmpath ../
end

function gam_o=SearchCritGamma(alp,bet,gam_min,gam_max)
tol=1e-6;
gam1=gam_min; 
[~,U1,~]=GetResponse(alp,bet,gam1,0);
isunder1=min(U1)<0; iter1=0; maxiter1=6;
while ~isunder1 && iter1<maxiter1
    gam1 = gam1/10;
    [~,U1,~]=GetResponse(alp,bet,gam1,0);
    isunder1=min(U1)<0;
    iter1=iter1+1;
end

gam2=gam_max;
[~,U2,~]=GetResponse(alp,bet,gam2,0);
isunder2=min(U2)<0; iter2=0; maxiter2=6;
while isunder2 && iter2<maxiter2
    gam2 = gam2*10;
    [~,U2,~]=GetResponse(alp,bet,gam2,0);
    isunder2=min(U2)<0;
    iter2=iter2+1;
end

if iter1<maxiter1  && iter2<maxiter2
    iter_c=0; maxiter_c=20;
    gamc = (gam1+gam2)/2;
    while isunder1 && ~isunder2 && abs(gam1-gam2)>tol && iter_c<maxiter_c
        [~,Uc,~]=GetResponse(alp,bet,gamc,0);
        isunderc=min(Uc)<0;
        if isunderc
            gam1=gamc;
            isunder1=isunderc;
        else
            gam2=gamc; 
            isunder2=isunderc;
        end
        gamc = (gam1+gam2)/2;
        iter_c=iter_c+1;
    end
    gam_o=gamc;
else
    if iter1>=maxiter1
        gam_0=gam1;
    elseif iter2>=maxiter2 
        gam_0=gam2;
    end
end
end

function [Tstore,Ustore,dUstore]=GetResponse(alp,bet,gam,gtil)
tspan=[0,100]; x0=[0;1]; Tstore=[]; Ustore=[]; Wstore=[];
Atol=1e-6; Rtol=1e-3; 
options=odeset('AbsTol',Atol,'RelTol',Rtol);
[T,X] = ode45(@(t,x)cor_ode_fun(t,x,gam,bet,alp,gtil),tspan,x0,options);
Tstore=[Tstore,T']; Ustore=[Ustore,X(:,1)']; Wstore=[Wstore,X(:,2)'];
m=length(Tstore);
dUstore=zeros(1,m); 
for i=1:length(Tstore)
    uplus=max([Ustore(i),0]);
    dUstore(i) = Wstore(i) - gam*uplus^bet;
end
end

function dx=cor_ode_fun(t,x,gam,bet,alp,gtil)
uplus = max([x(1,1),0]);
dx(1,1) = x(2,1) - gam*(uplus^bet);
dx(2,1) = - uplus^alp + gtil;
end
function gam_0=findcritgamma_0g(alp,bet,gam_a,gam_b)
tol=1e-3; gam_0=gam_b;
Ea=numericCOR(alp,bet,gam_a,0);
Eb=numericCOR(alp,bet,gam_b,0);
while abs(Eb - Ea)>=tol
   gam_0 = gam_b - Eb*(gam_b-gam_a)/(Eb - Ea);
   E0 = numericCOR(alp,bet,gam_0,0);
   gam_a=gam_b; Ea=Eb; 
   gam_b=gam_0; Eb=E0;
end
end