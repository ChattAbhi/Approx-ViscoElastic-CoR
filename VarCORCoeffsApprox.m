function varargout=VarCORCoeffsApprox(alp,bet,the)
% VarCORCoeffsApprox(alp,bet,the) computes the analytical approximations of
% $\mathcal{I}(\tilde{g}(\theta))$ and $\mathcal{I}(\tilde{g}(\theta))$, 
% when provided with the parameter values of alp, bet, and the, which 
% represent the parameters $\alpha$, $\beta$, and $\theta$, respectively.

[K,L]=computeKL(alp,bet,the);
J=Japprox(alp,bet,the);
M=Mapprox(alp,bet,the);
I=K*J; Q=L*M;

Jcubic=Japprox(alp,bet,the,'cubic');
Mfrac=Mapprox(alp,bet,the,'frac');
Mxlnx=Mapprox(alp,bet,the,'xlnx');
Mx2lnx=Mapprox(alp,bet,the,'x2lnx');


Icubic=K*Jcubic;
Qfrac=L*Mfrac; Qxlnx=L*Mxlnx; Qx2lnx=L*Mx2lnx; 

varargout{1} = I; 
varargout{2} = Q;
if nargout>=3 
    varargout{3} = Qfrac;
end
if nargout>=4
    varargout{4} = Qxlnx;
end
if nargout>=5
    varargout{5} = Qx2lnx; 
end
if nargout>=6
    varargout{6} = Icubic; 
end
    

end