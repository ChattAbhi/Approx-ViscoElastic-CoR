function [K,L]=computeKL(alp,bet,the)
% computeKL(alp,bet,the) computes the functions $K(\theta)$ and $L(\theta)$
% with the inputs alp, bet, and the, which stand for the parameters $\alpha$
% $\beta$, and $\theta$.
K=(2/sqrt(the))*(((alp+1)/(2*the))^(bet/(alp+1)));
L=sqrt(the)*((alp+1)/(2*the))^(bet/(alp+1));
end