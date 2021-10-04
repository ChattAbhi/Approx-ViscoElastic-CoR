function [I,Q]=VarCORCoeffsNumeric(alp,bet,gtil,varargin)
% VarCORCoeffsNumeric(alp,bet,gtil,_) computes the values of 
% $\mathcal{I}(\tilde{g})$ and $\matcal{Q}(\tilde{g})$ via numerical
% integration, for the given inputs of alp, bet, and gtil, which represent
% the parameters $\alpha$, $\beta$, and $\tilde{g}$, respectively.
% Additional arguments can be passed to set the Absolute and Relative Error
% Tolerences and Maximum Iterations for the numerical algorithm, as
% Name-Value pairs, with the names 'AbsTol', 'RelTol', and 'MaxIter',
% respectively.
y0=[0;0];
max_iter=100; uin=0; ufin=1; Atol=1e-8; Rtol=1e-7; Ustore=[]; Ystore=[];  
if ~isempty(varargin)
    for i=1:length(varargin)/2
        if strcmpi(varargin{2*i-1},'Abstol')
            Atol=varargin{2*i};
        end
        if strcmpi(varargin{2*i-1},'RelTol')
            Rtol=varargin{2*i};
        end
        if strcmpi(varargin{2*i-1},'MaxIter')
            max_iter=varargin{2*i};
        end
    end
end
iter=1;
while iter<=max_iter
    uspan=[uin,ufin];
    options=odeset('AbsTol',Atol,'RelTol',Rtol,'Event',@(u,y)coeff_event(u,y,alp,bet,gtil));
    [U,Y,Ue,Ye,Ie] = ode45(@(u,y)coeff_derivative(u,y,alp,bet,gtil),uspan,y0,options);
    uin=U(end); ufin=U(end)+1; y0=Y(end,:)'; Ustore=[Ustore,U']; Ystore=[Ystore,Y'];
    if ~isempty(Ie)
        break
    end
    iter=iter+1;
end
if iter>=max_iter
    error('Too many iterations increase Abstol and RelTol')
end
I=2*Ystore(1,end); Q=real(Ystore(2,end));
end
function dy=coeff_derivative(u,y,alp,bet,gtil)
% Derivative for numerical integration computation required for calculating
% I(gtilde) and Q(gtilde) the gtilde dependent coefficients for COR
dy=zeros(2,1);
dy(1) = (abs(1 + 2*(gtil*u - ((u^(alp+1))/(alp+ 1))))^(1/2))*(u^(bet-1)); %dI/du0
dy(2) = (u^(bet-1))*(abs(1 + 2*(gtil*u - ((u^(alp+1))/(alp+ 1))))^(-1/2)); %dQ/du0
end
function [val,ister,dir]=coeff_event(u,y,alp,bet,gtil)
% Event function for numerical integration of I(gtilde) and Q(gtilde) to stop 
% at maximum compression (The integrals are symmetric about max. compression 
% and therefore are reflected).
val=(1/(alp+1))*u^(alp+1) - gtil*u - 1/2;
ister=1; 
dir=0;
end