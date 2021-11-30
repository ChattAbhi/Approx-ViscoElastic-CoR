function gtil=gtilde_theta(alp,the)
% gtilde_theta(alp,the) computes the value of $\tilde{g}$ for a given
% values of $\theta$. The input alp and the represent the parameters 
% $\alpha$ and $\theta$ based on the viscoelastic model.

% Authors: Abhishek Chatterjee, Guillaume James, and Bernard Brogliato
% Address: Univ. Grenoble Alpes, INRIA, CNRS, Grenoble INP, LJK, Grenoble
%          38000 France 
gtil=(((alp + 1)/(2*the))^(-1/(alp+1)))*(1/(2*the) - 1/2);
end 