function Compare_Integ_Approx
prefixA='Figures/matfig/Components/integrand/'; prefixB='Figures/pdf/Components/integrand/';
the=[0,0.4,0.7,1];
alp=3/2;
N = [7,11,21,31];
x=0:0.01:1;
h=zeros(length(N)+1,length(x)); eh=zeros(length(N),length(x));
legA={}; legB={};
ls={'--r','-.g','--c','-.m','--k','-.y'};
for i=1:length(the)
    figA=figure('color','w'); figB=figure('color','w');
    AxA=axes(figA); AxB=axes(figB);
    for ii=1:length(x)
        h(1,ii)=func(alp,the(i),x(ii));
    end
    plot(AxA,x,h(1,:),'-b','LineWidth',1)
    legA{1}='$h_{\alpha,\theta}(x)$';
    for jj=1:length(N)
        n=N(jj);
        %afuns=lagrange_coeff_fun(n);
        c=ChebCoeff(n,alp,the(i));
        %a=afuns(alp,the(i));
        
        for ii=1:length(x)
            h(jj+1,ii) = 0; 
            for j=0:n
                h(jj+1,ii) = h(jj+1,ii) + c(j+1)*ChebT(j,x(ii)-1);
            end
%             h(2*jj,ii) = 0; h(2*jj+1,ii) = 0;
%             for j=0:n
%                 h(2*jj,ii) = h(2*jj,ii) + c(j+1)*ChebT(j,x(ii)-1);
%             end
%             for j=0:length(a)-1
%                 h(2*jj+1,ii) = h(2*jj+1,ii) + a(j+1)*x(ii)^j;
%             end
        end
        eh(jj,:)=abs(h(jj+1,:)-h(1,:));
        legA{jj+1}=strcat('$p_n(x)$, $n=',num2str(N(jj)),'$');
        legB{jj}=strcat('$p_n(x)$, $n=',num2str(N(jj)),'$');
        
%         eh(2*jj-1,:)=abs(h(2*jj,:)-h(1,:));
%         eh(2*jj,:)=abs(h(2*jj+1,:)-h(1,:));
%         legA{2*jj}=strcat('$p_n(x)$ - Chebyshev Basis, $n=',num2str(N(jj)),'$');
%         legA{2*jj+1}=strcat('$p_n(x)$ - Symbolic Lagrange, $n=',num2str(N(jj)),'$');
%         legB{2*jj-1}=strcat('$p_n(x)$ - Chebyshev Basis, $n=',num2str(N(jj)),'$');
%         legB{2*jj}=strcat('$p_n(x)$ - Symbolic Lagrange, $n=',num2str(N(jj)),'$');
        
        hold(AxA,'on') 
        plot(AxA,x,h(jj+1,:),ls{jj},'LineWidth',1) 
        hold(AxA,'on') 
        
        semilogy(AxB,x,eh(jj,:),ls{jj},'LineWidth',1)
        hold(AxB,'on')
        
%         hold(AxA,'on') 
%         plot(AxA,x,h(2*jj,:),ls{2*jj-1},'LineWidth',1) 
%         hold(AxA,'on')
%         plot(AxA,x,h(2*jj+1,:),ls{2*jj},'LineWidth',1)
%         hold(AxA,'on')
%         
%         semilogy(AxB,x,eh(2*jj-1,:),ls{2*jj-1},'LineWidth',1)
%         hold(AxB,'on')
%         semilogy(AxB,x,eh(2*jj,:),ls{2*jj},'LineWidth',1)
%         hold(AxB,'on')
    end
    
    hold(AxA,'off')
    legend(AxA,legA{:},'Interpreter','latex','FontSize',12,'Location','Best')
    ylabel(AxA,'$h_{\alpha,\theta}(x)$','Interpreter','latex','FontSize',16)
    xlabel(AxA,'$x$','Interpreter','latex','FontSize',16)
    title(AxA,strcat('$h_{\alpha,\theta}(x)$ \& $p_n(x)$ comparison for $\theta = ',num2str(the(i)),'$'),'Interpreter','latex','FontSize',16)
    savefig(figA,strcat(prefixA,'h_the',num2str(i),'.fig'));
    print(figA,strcat(prefixB,'h_the',num2str(i),'.pdf'),'-dpdf');
    
    hold(AxB,'off')
    legend(AxB,legB{:},'Interpreter','latex','FontSize',12,'Location','Best')
    ylabel(AxB,'Absolute Error, $h_{\alpha,\theta}(x)$','Interpreter','latex','FontSize',16)
    xlabel(AxB,'$x$','Interpreter','latex','FontSize',16)
    title(AxB,strcat('Absolute Error, $h_{\alpha,\theta}(x)$ \& $p_n(x)$ comparison for $\theta = ',num2str(the(i)), '$'),'Interpreter','latex','FontSize',16)
    savefig(figB,strcat(prefixA,'er_h_the',num2str(i),'.fig'));
    print(figB,strcat(prefixB,'er_h_the',num2str(i),'.pdf'),'-dpdf');
end



% legend('$h_{\alpha,\theta}(x)$','$p_n(x)$ - Chebyshev Basis, $n$=15', '$p_n(x)$ - Symbolic Lagrange, $n$=15','Interpreter','latex','FontSize',12)
% ylabel('$h_{\alpha,\theta}(x)$','Interpreter','latex','FontSize',16)
% xlabel('$x$','Interpreter','latex','FontSize',16)
% title(strcat('$h_{\alpha,\theta}(x)$ \& $p_n(x)$ comparison for $\theta = ',num2str(the(i)),'$ and $n = ',num2str(n),'$'),'Interpreter','latex','FontSize',16)
% 
% legend('$p_n(x)$ - Chebyshev Basis, $n$=15', '$p_n(x)$ - Symbolic Lagrange, $n$=15','Interpreter','latex','FontSize',12)
% ylabel('Absolute Error, $h_{\alpha,\theta}(x)$','Interpreter','latex','FontSize',16)
% xlabel('$x$','Interpreter','latex','FontSize',16)
% title(strcat('Absolute Error, $h_{\alpha,\theta}(x)$ \& $p_n(x)$ comparison for $\theta = ',num2str(the(i)), '$ and $n = ',num2str(n),'$'),'Interpreter','latex','FontSize',16)
end

function c=ChebCoeff(n,alp,the)
c = zeros(1,n+1);
for j=0:n
    c(j+1) = 0;
    for k = (n+1)/2:n
        xk = cos(pi*(2*k+1)/(2*n+2)) + 1;
        if j==0
            c(j+1) = c(j+1) + (2/(n+1))*func(alp,the,xk);
        elseif j>=2 && mod(j,2)==0
            c(j+1) = c(j+1) + (4/(n+1))*func(alp,the,xk)*cos(j*pi*(2*k+1)/(2*n+2));
        end
    end
end
end

function T=ChebT(j,x)
if j==0
    T=1;
elseif j==1 
    T=x;
elseif j>1
    T = 2*x*ChebT(j-1,x) - ChebT(j-2,x);
else
    error('The order of chebyshev polynomial must be a positive integer.')
end
end

function afuns=lagrange_coeff_fun(n)
alp=sym('alp','real'); 
the=sym('the','real');
x=sym('x','real');
len=length((n+1)/2:n);
xk=zeros(1,n+1);
yk=sym(zeros(1,n+1));
for k=0:n
    xk(k+1) = 1 + cos(pi*(2*k+1)/(2*n+2));
    if xk(k+1)<1
        yk(k+1) = func(alp,the,xk(k+1));
    else 
        yk(k+1) = func(alp,the,2-xk(k+1));
    end
end
af=sym(0);
for k=0:n
    xi = x*ones(1,n); xj=xk(k+1)*ones(1,n);
    nk=setdiff(1:n+1,k+1); 
    xm=xk(nk);
    af  = af + yk(k+1)*prod((xi-xm)./(xj-xm));
end

a=coeffs(af,x);
afuns=matlabFunction(a,'Vars',{alp,the});
end

function f=func(alp,the,x)
f = sqrt(((1-x^(alp+the))*(x^(1-the)))/(the + (1-the)*x - x^(alp+1)));
end