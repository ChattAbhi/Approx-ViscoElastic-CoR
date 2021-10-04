function varargout=numericCOR(alp,bet,gam,gtil,varargin)
% numericCOR(alp,bet,gam,gtil,_) numerically computes the Coefficient of
% Restitution, $e$. The input arguments alp, bet, gam, and gtil, represent
% the parameters $\alpha$, $\beta$, $gamma$, and $\tilde{g}$, repsectively.
% The Absolute and Relative Error Tolerences, and the Maximum Iteration
% Limit for the numerical algorithm can be specified as Name-Value pair
% argurments, using the names 'AbsTol', 'RelTol', and 'MaxIter'.
% Additionally, the evolution of the displacement $u(\tau)$ and the speed
% $u'(\tau)$ can be obtained by passing another Name-Value argument pairs
% using the name 'ShowPlots' as 
% numericCOR(alp,bet,gam,gtil,...,'ShowPlots',true). Save the plots as .fig
% and .pdf files by specifying the filenames using the names 'figsave' and
% 'pdfsave', respectively as, 
% numericCOR(...,'ShowPlots',true,'figsave','filename.fig',...
% 'pdfsave','filename.pdf'). The files would be saved as
% 'filename_resp.fig', 'filename_phase.fig', filename_resp.pdf', and 
% 'filename_phase.pdf'.

ShowPlots=false;
%Atol=1e-8; Rtol=1e-7; miter=100000; pdfname=[]; figname=[];
Atol=1e-14; Rtol=1e-10; miter=100000; pdfname=[]; figname=[]; mult=false; ignore=false; outDat=false;
if ~isempty(varargin)
    for i=1:length(varargin)/2
        if strcmpi(varargin{2*i-1},'Abstol')
            Atol=varargin{2*i};
        end
        if strcmpi(varargin{2*i-1},'RelTol')
            Rtol=varargin{2*i};
        end
        if strcmpi(varargin{2*i-1},'MaxIter')
            miter=varargin{2*i};
        end
        if strcmpi(varargin{2*i-1},'MultiplyIter')
            mult=varargin{2*i};
        end
        if strcmpi(varargin{2*i-1},'IgnoreMaxIterError')
            ignore=varargin{2*i};
        end
        if strcmpi(varargin{2*i-1},'ShowPlots')
            ShowPlots=varargin{2*i};
        end
        if strcmpi(varargin{2*i-1},'pdfsave')
            pdfname=varargin{2*i};
        end
        if strcmpi(varargin{2*i-1},'figsave')
            figname=varargin{2*i};
        end
        if strcmpi(varargin{2*i-1},'OutData')
            outDat=varargin{2*i};
        end
    end
end

max_iter=miter; iter=0; tin=0; tfin=100; x0=[0;1]; Tstore=[]; Ustore=[]; Wstore=[];
inter=1;
while iter<=max_iter
    tspan=[tin,tfin];
    options=odeset('AbsTol',Atol,'RelTol',Rtol,'Event',@(t,x)cor_event(t,x,gam,bet,alp,gtil,inter));
    [T,X,Te,Xe,Ie] = ode45(@(t,x)cor_ode_fun(t,x,gam,bet,alp,gtil,inter),tspan,x0,options);
    tin=T(end); tfin=T(end)+100; x0=X(end,:)'; Tstore=[Tstore,T']; Ustore=[Ustore,X(:,1)']; Wstore=[Wstore,X(:,2)'];
    if mult
        tfin=max([T(end)+10,100*T(end),T(end)*exp(0.1*iter)]);
    end
    if ~isempty(Ie)
        inter = inter + 1;
    end
    if inter==3 && (X(end,1)<=0 || X(end,2) - gam*(max([X(end,1),0]))^bet>=0)
        break;
    end
    if inter>3
        break;
    end
    iter=iter+1;
end
if iter>=max_iter && ~ignore
    error('Too many iterations increase Abstol and/or RelTol!')
end
m=length(Tstore);
dUstore=zeros(1,m); 
for i=1:length(Tstore)
    uplus=max([Ustore(i),0]);
    dUstore(i) = Wstore(i) - gam*uplus^bet;
end

e = -dUstore(end);

if outDat
    varargout{1}=e; varargout{2}=Tstore; varargout{3}=Ustore; varargout{4}=dUstore;
else
    varargout{1}=e;
end

if ShowPlots
    figure('Color','w')
    hold off
    plot([Tstore(1),1.1*Tstore(end)],[Ustore(end),Ustore(end)],'-k')
    hold on
    plot([Tstore(1),1.1*Tstore(end)],[-e,-e],'-k')
    hold on
    plt_u=plot(Tstore,Ustore,'-b','LineWidth',1);
    hold on
    plt_up=plot(Tstore,dUstore,'--r','Linewidth',1);
    xlim([Tstore(1),1.1*Tstore(end)])
    xlabel('$$\tau$$','FontSize',16,'Interpreter','Latex')
    ylabel('$$u(\tau)$$,$$u^{\prime}(\tau)$$','FontSize',16,'Interpreter','latex')
    legend([plt_u,plt_up],{'$$u(\tau)$$','$$u^{\prime}(\tau)$$'},'FontSize',16,'Interpreter','latex','Location','Best')
    title(strcat('$$\alpha = ',num2str(alp,2),'$$, $$\beta = ',num2str(bet,2),'$$, $$\gamma = ',num2str(gam,4),'$$, and $$\tilde{g} = ',num2str(gtil,3),'$$'),'Interpreter','latex','FontSize',16)
    text(3/4*Tstore(end),Ustore(end),strcat('$u=',num2str(round(Ustore(end),4,'decimals')),'$'),'Interpreter','latex','FontSize',16)
    text(Tstore(end)*(0.95),-0.95*e,strcat('$$e=',num2str(round(e,4,'decimals')),'$$'),'Interpreter','latex','FontSize',16)
    if ~isempty(figname)
        name_comp=strsplit(figname,'.');
        name_parts={name_comp{1},'_resp.',name_comp{2}};
        newfigname=strcat(name_parts{:});
        savefig(newfigname);
    end
    if ~isempty(pdfname)
        name_comp=strsplit(pdfname,'.');
        name_parts={name_comp{1},'_resp.',name_comp{2}};
        newpdfname=strcat(name_parts{:});
        print(newpdfname,'-dpdf');
    end
    
    figure('Color','w')
    plot([0,0],[1.1*min(dUstore),1.1*max(dUstore)],'-k')
    hold on 
    plot([1.1*min(Ustore),1.1*max(Ustore)],[0,0],'-k')
    hold on
    plot([0,0],[1.1*min(dUstore),0],'-r','LineWidth',2)
    hold on
    plot([0,gtil^(1/alp)],[0,0],'-r','LineWidth',2)
    hold on
    plot(gtil^(1/alp),0,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5)
    hold on
    plot(Ustore,dUstore,'-b','LineWidth',1)
    axis([1.1*min(Ustore),1.1*max(Ustore),1.1*min(dUstore),1.1*max(dUstore)])
    xlabel('$$u(\tau)$$','FontSize',16,'Interpreter','Latex')
    ylabel('$$u^{\prime}(\tau)$$','FontSize',16,'Interpreter','latex')
    title(strcat('$$\alpha = ',num2str(alp,2),'$$, $$\beta = ',num2str(bet,2),'$$, $$\gamma = ',num2str(gam,4),'$$, and $$\tilde{g} = ',num2str(gtil,3),'$$'),'Interpreter','latex','FontSize',16)
    text(gtil^(1/alp),0.1*max(dUstore),strcat('$\bar{u}={\tilde{g}}^{\frac{1}{\alpha}}=',num2str(round(gtil^(1/alp),4,'decimal')),'$'),'Interpreter','latex','FontSize',16)
    grid on
    if ~isempty(figname)
        name_comp=strsplit(figname,'.');
        name_parts={name_comp{1},'_phase.',name_comp{2}};
        newfigname=strcat(name_parts{:});
        savefig(newfigname);
    end
    if ~isempty(pdfname)
        name_comp=strsplit(pdfname,'.');
        name_parts={name_comp{1},'_phase.',name_comp{2}};
        newpdfname=strcat(name_parts{:});
        print(newpdfname,'-dpdf');
    end
    
    fprintf('e = %.4f\n',e)
    fprintf('g^(1/alpha) = %.4f\n',gtil^(1/alp))
    fprintf('u_end = %.4f\n',Ustore(end))
end
end

function dx=cor_ode_fun(t,x,gam,bet,alp,gtil,inter)
uplus = max([x(1,1),0]);
dx(1,1) = x(2,1) - gam*(uplus^bet);
dx(2,1) = - uplus^alp + gtil;
end

function [val,ister,dir]=cor_event(t,x,gam,beta,alp,gtil,inter)
u=x(1,1);
uplus=max([u,0]);
udot=x(2) - gam*uplus^beta;
if inter==1
    val=[udot;1;1];
    ister=[1;0;0];
    dir=[-1;0;0];
elseif inter==2 
    val=[u-gtil^(1/alp);1;1];
    ister=[1;0;0];
    dir=[-1;0;0];
else
    val=[u*udot;u;udot];
    ister=[1;1;1];
    dir=[1;-1;1];
end
end


