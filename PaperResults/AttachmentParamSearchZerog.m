function AttachmentParamSearchZerog()
addpath ../

alp1=1; alp2=5/4; alp3=3/2; alp4=5/2; bet=1:0.05:3;
gam_c=zeros(4,length(bet));
for i=1:length(bet)
    %gam_c(i)=findzero_e_gam(alp,bet(i),1e-4,0.1);
    f1=findcritgamma_0g(alp1,bet(i),1e-4,0.1);
    %f1=findzero_e_gam(alp1,bet(i),1e-4,0.1);
    if f1>10
        gam_c(1,i)=nan;
    else
        gam_c(1,i)=f1;
    end
    f2=findcritgamma_0g(alp2,bet(i),1e-4,0.1);
    %f2=findzero_e_gam(alp2,bet(i),1e-4,0.1);
    if f1>10
        gam_c(2,i)=nan;
    else
        gam_c(2,i)=f2;
    end
    f3=findcritgamma_0g(alp3,bet(i),1e-4,0.1);
    %f3=findzero_e_gam(alp3,bet(i),1e-4,0.1);
    if f3>10
        gam_c(3,i)=nan;
    else
        gam_c(3,i)=f3;
    end
    f4=findcritgamma_0g(alp4,bet(i),1e-4,0.1);
    %f4=findzero_e_gam(alp4,bet(i),1e-4,0.1);
    if f4>10
        gam_c(4,i)=nan;
    else
        gam_c(4,i)=f4;
    end
%     gam_c(1,i)=findcritgamma_0g(alp1,bet(i),1e-4,0.1);
%     gam_c(2,i)=max([10,findcritgamma_0g(alp2,bet(i),1e-4,0.1)]);
%     gam_c(3,i)=max([10,findcritgamma_0g(alp3,bet(i),1e-4,0.1)]);
%     gam_c(4,i)=findcritgamma_0g(alp4,bet(i),1e-4,0.1);
    fprintf('beta = %.4f\n',bet(i))
end
figure();
plot(bet,gam_c(1,:),'-b',bet,gam_c(2,:),'-r',bet,gam_c(3,:),'-g',bet,gam_c(4,:),'-m')
legend('$\alpha=1$','$\alpha=5/4$','$\alpha=3/2$','$\alpha=5/2$','FontSize',20,'Interpreter','latex','Location','Best')
xlabel('$\beta$','FontSize',20,'Interpreter','latex')
ylabel('$\gamma_c$','FontSize',12,'Interpreter','latex')

figure();
semilogy(bet,gam_c(1,:),'-b',bet,gam_c(2,:),'-r',bet,gam_c(3,:),'-g',bet,gam_c(4,:),'-m')
legend('$\alpha=1$','$\alpha=5/4$','$\alpha=3/2$','$\alpha=5/2$','FontSize',20,'Interpreter','latex','Location','Best')
xlabel('$\beta$','FontSize',20,'Interpreter','latex')
ylabel('$\gamma_c$','FontSize',20,'Interpreter','latex')



alp=1:0.05:3; bet=1:0.05:3;
%alp=logspace(1,2,100); bet=logspace(1,2,100);
%gam_c=zeros(length(alp),length(bet));
[Alp,Bet]=meshgrid(alp,bet); 
[m,n]=size(Alp);
gam_c=zeros(m,n);
for i=1:m
    for j=1:n
        %gam_c(i,j)=findzero_e_gam(Alp(i,j),Bet(i,j),1e-4,0.1);
        gam_c(i,j)=findcritgamma_0g(Alp(i,j),Bet(i,j),1e-4,0.1);
        fprintf('[i,j]=[%i,%i]\n',i,j)
    end
end

figure('Color','w')
indx1=find(alp==1); indx2=find(alp==5/4); indx3=find(alp==3/2); indx4=find(alp==5/2);
semilogy(Bet(:,indx1),gam_c(:,indx1),'-b')
hold on 
semilogy(Bet(:,indx2),gam_c(:,indx2),'-r')
hold on 
semilogy(Bet(:,indx3),gam_c(:,indx3),'-g')
hold on 
semilogy(Bet(:,indx4),gam_c(:,indx4),'-m')
legend('$\alpha=1$','$\alpha=5/4$','$\alpha=3/2$','$\alpha=5/2$','FontSize',20,'Interpreter','latex','Location','Best')
xlabel('$\beta$','FontSize',20,'Interpreter','latex')
ylabel('$\gamma_c$','FontSize',20,'Interpreter','latex')
title(strcat('$\gamma_c$ vs. $\beta$, with $\tilde{g} = 0$'), 'FontSize', 20, 'Interpreter', 'latex')
savefig('Figures/matfig/NoReboundTerm/zeroG_no_term_param.fig');
print('Figures/pdf/NoReboundTerm/zeroG_no_term_param.pdf', '-dpdf');

figure('Color','w')
surf(Alp,Bet,log10(gam_c))
colorbar
xlabel('$\alpha$','FontSize',20,'Interpreter','latex')
ylabel('$\beta$','FontSize',20,'Interpreter','latex')
zlabel('$\log_{10}  \gamma_c $','FontSize',20,'Interpreter','latex')
title('$\log_{10} \gamma_c (\alpha,\beta)$, when $\tilde{g}=0$','FontSize',20,'Interpreter','latex')
savefig('Figures/matfig/NoReboundTerm/log_surf_zeroG_no_term_param.fig');
print('Figures/pdf/NoReboundTerm/log_surf_zeroG_no_term_param.pdf', '-dpdf');

figure('Color','w')
contour(Alp,Bet,log10(gam_c),'ShowText','on')
colorbar
xlabel('$\alpha$','FontSize',20,'Interpreter','latex')
ylabel('$\beta$','FontSize',20,'Interpreter','latex')
title('$\log_{10} \gamma_c (\alpha,\beta)$, when $\tilde{g}=0$','FontSize',20,'Interpreter','latex')

figure('Color','w')
contourf(Alp,Bet,log10(gam_c),'ShowText','on')
colorbar
xlabel('$\alpha$','FontSize',20,'Interpreter','latex')
ylabel('$\beta$','FontSize',20,'Interpreter','latex')
title('$\log_{10} \gamma_c (\alpha,\beta)$, when $\tilde{g}=0$','FontSize',20,'Interpreter','latex')
savefig('Figures/matfig/NoReboundTerm/log_cont_zeroG_no_term_param.fig');
print('Figures/pdf/NoReboundTerm/log_cont_zeroG_no_term_param.pdf', '-dpdf');

%numericCOR(3/2,5/2,gam_0,0)
rmpath ../
end

function gam_0=findzero_e_gam(alp,bet,gam_min,gam_max)
max_iter=100;
iter=0; 

ztol=1e-4; rtol=1e-3; 

gam_a=gam_min; gam_b=gam_max; range_search = true;
gam_c=gam_b;
while iter<=max_iter
    Eb=numericCOR(alp,bet,gam_b,0,'AbsTol',1e-6,'RelTol',1e-3,'MaxIter',10000,'MultiplyIter',true);
    if range_search
        if Eb<ztol
            range_search=false;           
        else
            %gam_b=gam_b*((iter+1)^2);
            gam_b=2*gam_b*exp(0.01*iter);
        end
    else
        gam_c = (gam_a + gam_b)/2;
        Ec = numericCOR(alp,bet,gam_c,0,'AbsTol',1e-6,'RelTol',1e-3,'MaxIter',10000,'MultiplyIter',true);
        if Ec<ztol
            gam_b=gam_c;
        else
            gam_a=gam_c;
        end
    end 
    if abs(gam_a-gam_b)<=rtol
        gam_0 = (gam_a + gam_b)/2;
        E0 = numericCOR(alp,bet,gam_0,0,'AbsTol',1e-6,'RelTol',1e-3,'MaxIter',10000,'MultiplyIter',true);
        if E0>=ztol
            gam_0=gam_b;
        end
        break;
    end
    iter=iter+1;
end
if iter>max_iter
    %error('Too many iterations. No solution found.')
    gam_0=nan;
end
close all
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