function J=Japprox(alp,bet,the,varargin)
% Japprox(alp,bet,the) computes the analytical (cubic-Hermite) approximation 
% of $J_{\alpha,\beta}(\theta)$, with the inputs of alp, bet and the, which 
% respectively repersent the parameters $\alpha$, $\beta$, and $\theta$.

sel='M'; 
if ~isempty(varargin)
    sel=varargin{1};
end

if strcmpi(sel,'cubic')
    p1=(1/(alp+1))*beta(bet/(alp+1),3/2);
    m1=(1/(2*(alp+1)))*(beta(bet/(alp+1),1/2) - beta((bet + 1)/(alp + 1),1/2));
    p0=(1/alp)*beta((2*bet + 1)/(2*alp),3/2);
    m0=(1/(2*alp))*( beta((2*bet-1)/(2*alp),1/2) - beta((2*bet+1)/(2*alp),1/2) );
    
    J=(2*the^3 - 3*the^2 + 1)*p0 + (the^3 - 2*the^2 + the)*m0 + (-2*the^3 + 3*the^2)*p1 + (the^3 - the^2)*m1;
elseif strcmpi(sel,'M')
    Ma=Mapprox(alp,alp+bet+1,the);
    Mb=Mapprox(alp,bet+1,the);
    
    J = ((alp+1)/(2*bet))*Ma + ((the-1)/(2*bet))*Mb;
else
    error('!!')
end
end