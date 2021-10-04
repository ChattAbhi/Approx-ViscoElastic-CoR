function gamc=gamma_c(alp,bet,gtil)
% gam_c(alp,bet,gtil) returns the value of gamma such that $e=0$, with the
% arguments alp, bet, and gtil, which represent the parmeters $\alpha$, 
% $\beta$, and $\tilde{g}$, respectively.
C=2*sqrt(2)*(bet/alp)*((alp+1)^(bet/alp + 1/(2*alp)))*beta((bet + 1/2)/alp, 3/2);
r=bet/alp + 1/(2*alp) + 1/2;
gamc=1/(2*C*(gtil^r));
end