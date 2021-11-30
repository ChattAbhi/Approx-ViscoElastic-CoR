function J=Japprox(alp,bet,the,varargin)
% Japprox(alp,bet,the) computes the analytical  approximation of the 
% integral $J_{\alpha,\beta}(\theta)$, with the inputs of alp, bet and the, 
% which respectively represent the parameters $\alpha$, $\beta$, and 
% $\theta$, used in the viscoelastic model.

% Authors: Abhishek Chatterjee, Guillaume James, and Bernard Brogliato
% Address: Univ. Grenoble Alpes, INRIA, CNRS, Grenoble INP, LJK, Grenoble
%          38000 France 

Ma=Mapprox(alp,alp+bet+1,the);
Mb=Mapprox(alp,bet+1,the);
J = ((alp+1)/(2*bet))*Ma + ((the-1)/(2*bet))*Mb;
end