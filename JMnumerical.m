function [J,M]=JMnumerical(alp,bet,the,varargin)
% JMnumerical(alp,bet,the,_) numerically computes the integrals
% $J_{\alpha,\beta}(\theta)$ and $M_{\alpha,\beta}(\theta)$, using the input
% arguments alp, bet, and the, which represent the parameters $\alpha$, 
% $\beta$, and $\theta$. Additionally, the Absolute and Relative Error 
% Tolerences can be specified as Name-Value pairs with the names 'AbsTol' 
% and 'RelTol', respectively.
Atol=1e-6; Rtol=1e-3;
if ~isempty(varargin)
    for i=1:length(varargin)/2
        if strcmpi(varargin{2*i-1},'Abstol')
            Atol=varargin{2*i};
        end
        if strcmpi(varargin{2*i-1},'RelTol')
            Rtol=varargin{2*i};
        end
    end
end
M=integral(@(x)funcM(x,the,alp,bet),0,1);
J=integral(@(x)funcJ(x,the,alp,bet),0,1);
end

function dM=funcM(x,th,alp,bet)
dM=((th + (1-th).*x - x.^(alp+1)).^(-1/2)).*(x.^(bet-1));
end

function dJ=funcJ(x,th,alp,bet)
dJ=((th + (1-th).*x - x.^(alp+1)).^(1/2)).*(x.^(bet-1));
end