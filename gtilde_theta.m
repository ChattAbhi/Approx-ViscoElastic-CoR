function gtil=gtilde_theta(alp,the)
% gtilde_theta(alp,the) computes the value of $\tilde{g}$ for a given
% values of $\theta$. The inputs alp and the represent the parameters 
% $\alpha$ and $\beta$
gtil=(((alp + 1)/(2*the))^(-1/(alp+1)))*(1/(2*the) - 1/2);
end 