function e=analyticCOR(alp,bet,gam,gtil)
% analyticCOR(alp,bet,gam,gtil) computes the analytical Coefficient of
% Restitution result, $e$. The input arguments alp, bet, gam, and gtil,
% represent the parameters $\alpha$, $\beta$, $\gamma$, and $\tilde{g}$,
% respectively.
the=theta_gtilde(alp,gtil);
[I,Q]=VarCORCoeffsApprox(alp,bet,the);
gam_c=2/(bet*(Q-(1/2)*I));

if gam>gam_c
    esqr=1 - 2*bet*I*gam + 2*(bet^2)*I*Q*(gam^2);
    e = sqrt(max([0,esqr])); 
else
    e = 1 - bet*I*gam + (bet^2)*I*(Q - (1/2)*I)*(gam^2);
end

end