%This file presents a simulation example of a bouncing steel ball when
%dropped from an initial height (datum is considered at the bottom most 
%point of the ball). The coefficient of restituion for each impact is 
%computed using one of the proposed methods from the paper (different 
%methods from the paper can be selected by modifying the variable "method"
%in this code). This code provides a working example of how the getCOR()
%function in the library, can be used to compute the CoR of a bead during
%simulation.
%
% Authors: Abhishek Chatterjee, Guillaume James, and Bernard Brogliato
% Address: Univ. Grenoble Alpes, INRIA, CNRS, Grenoble INP, LJK, Grenoble
%          38000 France 

clear all; close all; clc;

addpath ../

m=1.54e-1; k=3.6138e10; gam0=1.5237e-6; % Steel ball model parameters
alp=3/2; bet=3/2; % Viscoelastic contact model parameters
h0=1; v0=0; % Initial height and velocity
g=9.8; % gravitational acceleration constant

method='ord2approx'; % Choice of method for computing the CoR. The value 
                     % of this variable can be changed to other valid 
                     % method options, which are 'ord1approx',
                     % 'ord2sch-pos'(valid only if bet=3/2), 'ord2numint', 
                     % and 'dirnumint'. Type >> help getCOR in the command 
                     % window for more details on each method.

                    
tin=0; tfin=10; % Initial and final times for the simulation 
thr=1e-6; ntsp=100; % Attachment velocity threshold and the number of 
                    % output time-steps between each impact

repeat=true;
t0=tin; 
tstore=[]; hstore=[]; vstore=[]; estore=[]; % Simulation data storage
while repeat
    tn = t0 + max([real((v0 + sqrt(v0^2 + 4*g*h0))/(2*g)),real((v0 - sqrt(v0^2 + 4*g*h0))/(2*g))]);
    if tn>=tfin
        tn=tfin;
        repeat=false;
    end
    hn = h0 + v0*(tn-t0) - g*(tn-t0)^2;
    vn = v0 - 2*g*(tn-t0);
    ts=linspace(t0,tn,ntsp);
    tsd=linspace(0,tn-t0,ntsp);
    hs=h0*ones(1,ntsp) + v0.*tsd - g.*(tsd.^2);
    vs=v0*ones(1,ntsp) - 2.*g.*tsd;
    tstore=[tstore,ts]; hstore=[hstore,hs]; vstore=[vstore,vs];
    if repeat
        if abs(vn)>thr
            e=getCOR(m,k,gam0,-vn,alp,bet,g,0,method); % Computes the CoR. 
            %Note that pre-impact velocity passed to getCoR() function
            %must have a positive value.
            v0=-e*vn; % Post-impact velocity using the computed CoR.
            t0=tn; h0=hn;
            estore=[estore,e];
        else
            v0=0;
            t0=tn;
            h0=hn;
            ts=linspace(t0,tfin,ntsp);
            hs=h0*ones(1,ntsp);
            vs=v0*ones(1,ntsp);
            tstore=[tstore,ts]; hstore=[hstore,hs]; vstore=[vstore,vs];
            repeat=false;
        end
    end
end

% Plot height w.r.t. time
figure('color','w')
plot(tstore,hstore,'-b','Linewidth',1)
ylabel('Height, $h(t)[m]$','Fontsize',16,'Interpreter','latex')
xlabel('Time, $t[s]$','Fontsize',16,'Interpreter','latex')
title(strcat('Height of a bouncing steel ball, given $\alpha = ',...
    num2str(alp),'$ and $\beta = ',num2str(bet),'$'),'Fontsize',16,'Interpreter','latex');
axis([tin,tfin,0,1])
print('bounce_height','-dpng')

% Plot velocity w.r.t. time
figure('color','w')
plot(tstore,vstore,'-b','Linewidth',1)
ylabel('Velocity, $v(t)[m/s]$','FontSize',16,'Interpreter','latex')
xlabel('Time, $t[s]$','FontSize',16,'Interpreter','latex')
title(strcat('Velocity of a bouncing steel ball, given $\alpha = ',...
    num2str(alp),'$ and $\beta = ',num2str(bet),'$'),'Fontsize',16,'Interpreter','latex');
axis([tin,tfin,min(vstore),max(vstore)])
print('bounce_velocity','-dpng')

% Bar plot of CoR for each impact
figure('color','w')
bar(1:length(estore),estore)
ylabel('CoR, $e$','Fontsize',16,'Interpreter','latex')
xlabel('Impact Number','Fontsize',16,'Interpreter','latex')
title('Coefficient of Restitution value for each impact','Fontsize',16,'Interpreter','latex')
print('bounce_CoR','-dpng')

rmpath ../