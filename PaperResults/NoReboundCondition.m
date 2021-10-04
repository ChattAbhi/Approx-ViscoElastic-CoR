function NoReboundCondition(alp,bet,varargin)
gam_min = 0.001;
gam_max = 0.1;
gam_max_z = 0.01;
fname1 = 'Stopg'; fname2 = 'Stopg_err';
if ~isempty(varargin)
    if ~isempty(varargin{1})
        gam_min = varargin{1};
    end
    if ~isempty(varargin{2})
        gam_max = varargin{2};
    end
    if ~isempty(varargin{3})
        gam_max_z = varargin{3};
    end
    if ~isempty(varargin{4})
        fname1 = varargin{4};
    end
    if ~isempty(varargin{5})
        fname2 = varargin{5};
    end
end

num_steps=100;

gtilde_min=0;
gtilde_max=10;

gam_range =  linspace(gam_min,gam_max,num_steps);
%gam_range =  logspace(gam_min,gam_max,num_steps);
gtnum = zeros(1,length(gam_range)); 
gtanl = zeros(1,length(gam_range));
relerr = zeros(1,length(gam_range));
for i=1:length(gam_range)
    try
        gtnum(i)=findminzero_e(alp,bet,gam_range(i),gtilde_min,gtilde_max);
    catch
        gtnum(i)=nan;
    end
    %gtnum(i)=findgtilc(alp,bet,gam_range(i),0,0.01);
    C=2*sqrt(2)*(bet/alp)*((alp+1)^(bet/alp + 1/(2*alp)))*beta((bet+1/2)/alp,3/2);
    B=2*gam_range(i)*C; A=(bet/alp + 1/(2*alp) + 1/2);
    gtanl(i)=(1/B)^(1/A);
    relerr(i)=abs(gtnum(i)-gtanl(i))/gtnum(i);
    fprintf('gamma = %.4f, gtilde_num = %.4f, gtilde_anl = %.4f \n',gam_range(i),gtnum(i),gtanl(i));
end

gam_range_z = linspace(gam_min,gam_max_z,num_steps);
gtnumz = zeros(1,length(gam_range_z));
gtanlz = zeros(1,length(gam_range_z));
relerrz = zeros(1,length(gam_range_z));
for i = 1:length(gam_range_z)
    gtnumz(i) = findminzero_e(alp,bet,gam_range_z(i),gtilde_min,gtilde_max);
    %gtnumz(i)=findgtilc(alp,bet,gam_range_z(i),0,0.01);
    C=2*sqrt(2)*(bet/alp)*((alp+1)^(bet/alp + 1/(2*alp)))*beta((bet+1/2)/alp,3/2);
    B=2*gam_range_z(i)*C; A=(bet/alp + 1/(2*alp) + 1/2);
    gtanlz(i) = (1/B)^(1/A);
    relerrz(i) = abs(gtnumz(i) - gtanlz(i))/gtnumz(i);
    fprintf('gamma = %.4f, gtilde_num = %.4f, gtilde_anl = %.4f \n',gam_range_z(i),gtnumz(i),gtanlz(i));
end

fig=figure('Color','w'); 
ax_base=axes(fig);
h(1) = plot(ax_base,gam_range,gtnum,'-b','LineWidth',1);
hold on
h(2) = plot(ax_base,gam_range,gtanl,'--r','LineWidth',1);
hold off
grid(ax_base,'on')
title(strcat('$\tilde{g}_c$ vs. $\gamma$, with $\alpha = ',num2str(alp),'$ \& $\beta = ',num2str(bet),'$'),'Interpreter','Latex','FontSize',20)
legend(h,{'Numerical $\tilde{g}_c$','Analytical $\tilde{g}_c$, [Eq.~24]'},'Interpreter','latex','FontSize',20)
xlabel('$\gamma$','Interpreter','latex','FontSize',20)
ylabel('$\tilde{g}_c$','Interpreter','latex','FontSize',20)

pos=get(ax_base,'Position');
scl=0.65; dims=[pos(3)-pos(1),pos(4)-pos(2)];
ax_zoom=axes(fig,'OuterPosition',[pos(1),pos(2),scl*dims(1),scl*dims(2)],'Box','on');

h(1) = plot(ax_zoom,gam_range_z,gtnumz,'-b','LineWidth',1);
hold on
h(2) = plot(ax_zoom,gam_range_z,gtanlz,'--r','LineWidth',1);
hold off
grid(ax_zoom,'on');

savefig(strcat('Figures/matfig/NoReboundTerm/',fname1,'.fig'));
%print(strcat('Figures/pdf/NoReboundTerm/',fname1,'.pdf'), '-dpdf');

fig=figure('Color','w');
ax_base=axes(fig);
plot(ax_base,gam_range(find(gam_range<=0.3)),relerr(find(gam_range<=0.3)),'-b','LineWidth',1);
grid(ax_base,'on')
title(strcat('Relative Error in $\tilde{g}_c$, with $\alpha = ',num2str(alp),'$ \& $\beta =',num2str(bet),'$'),'Interpreter','Latex','FontSize',20)
xlabel('$\gamma$','Interpreter','latex','FontSize',20)
ylabel('Relative Error','Interpreter','latex','FontSize',20)

pos=get(ax_base,'Position');
scl=0.65; dims=[pos(3)-pos(1),pos(4)-pos(2)];
ax_zoom=axes(fig,'OuterPosition',[pos(1),pos(2),scl*dims(1),scl*dims(2)],'Box','on');
plot(ax_zoom,gam_range_z,relerrz,'-b','LineWidth',1);
grid(ax_zoom,'on')

savefig(strcat('Figures/matfig/NoReboundTerm/',fname2,'.fig'));
%print(strcat('Figures/pdf/NoReboundTerm/',fname2,'.pdf'), '-dpdf');
end

function gt0=findminzero_e(alp,bet,gam,gtilde_min,gtilde_max)
max_iter=100;
iter=0;

ztol=1e-4; rtol=1e-3; 

gta=gtilde_min; gtb=gtilde_max; range_search = true;
gtc=gtb;
while iter<=max_iter
    %Ea=numericCOR(alp,bet,gam,gta);
    Eb=numericCOR(alp,bet,gam,gtb,'MaxIter',20,'IgnoreMaxIterError',true);
    if range_search
        if Eb<ztol
            range_search=false;           
        else
            gtb = gtb + 10;
        end
    else
        gtc = (gta + gtb)/2;
        Ec = numericCOR(alp,bet,gam,gtc,'MaxIter',20,'IgnoreMaxIterError',true);
        if Ec<ztol
            gtb=gtc;
        else
            gta=gtc;
        end
    end 
    if abs(gta-gtb)<=rtol
        gt0 = (gta + gtb)/2;
        E0 = numericCOR(alp,bet,gam,gt0,'MaxIter',20,'IgnoreMaxIterError',true);
        if E0>=ztol
            gt0=gtb;
        end
        break;
    end
    iter=iter+1;
end
if iter>max_iter
    error('Too many iterations. No solution found.')
end
end

function gtilc=findgtilc(alp,bet,gam,gtila,gtilb)
tol=1e-5; gtilc=gtilb;
Ea=numericCOR(alp,bet,gam,gtila);
Eb=numericCOR(alp,bet,gam,gtilb);
while abs(Eb-Ea)>=tol
    gtilc = gtilb - 0.1*Eb*(gtilb-gtila)/(Eb-Ea);
    E0 = numericCOR(alp,bet,gam,gtilc);
    gtila = gtilb; Ea=Eb;
    gtilb = gtilc; Eb=E0;
end
end