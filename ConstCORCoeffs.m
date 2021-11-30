function varargout=ConstCORCoeffs(alp,bet)
% ConstCORCoeffs(alp,bet) computes the constant coefficients in the
% expansion of $e=1 - \gamma C_0 - \gamma \tilde{g} C_1 - \gamma^2 C_2$.
% The input arguments alp and bet represent the parameters $\alpha$ and
% $\beta$. The output C2 which stand for the coefficient $C_2$ is only 
% computed for the valid case $\alpha = \beta = 3/2$.

% Authors: Abhishek Chatterjee, Guillaume James, and Bernard Brogliato
% Address: Univ. Grenoble Alpes, INRIA, CNRS, Grenoble INP, LJK, Grenoble
%          38000 France 

c=0.7982665553; %Constant from Schwager-Poschel
C0 = bet*(((alp+1)/2)^((bet-alp-1)/(alp+1)))*beta(bet/(alp+1),3/2);  % Coeff. for O(gtilde) term
C1 = bet*(((alp+1)/2)^((bet-alp)/(alp+1)))*beta((bet+1)/(alp+1),1/2); % Coeff. for O(gamma) term 
C2 = -c*(9/4); % Coeff. for O(gamma^2) term, found using Schwager-Poschel constant (valid for alp=bet=3/2)
if nargout==3 && bet==3/2
    varargout={C0,C1,C2};
elseif nargout==3 && bet~=3/2
    error('Second order constant coefficient is not known when beta~=3/2!')
elseif nargout<=2 && nargout>0
    varargout={C0,C1};
else
    error('Too many outputs!')
end
end