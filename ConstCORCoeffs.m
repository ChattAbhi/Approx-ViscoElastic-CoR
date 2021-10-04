function [C0,C1,C2]=ConstCORCoeffs(alp,bet)
% ConstCORCoeffs(alp,bet) computes the constant coefficients in the
% expansion of $e=1 - \gamma C_0 - \gamma \tilde{g} C_1 - \gamma^2 C_2$.
% The input arguments alp and bet represent the parameters $\alpha$ and
% $\beta$. The output C2 which stand for the coefficient $C_2$ is only 
% valid for the case $\alpha = \beta = 3/2$.

c=0.7982665553; %Constant from Schwager-Poschel
C0 = bet*(((alp+1)/2)^((bet-alp-1)/(alp+1)))*beta(bet/(alp+1),3/2);  % Coeff. for O(gtilde) term
C1 = bet*(((alp+1)/2)^((bet-alp)/(alp+1)))*beta((bet+1)/(alp+1),1/2); % Coeff. for O(gamma) term 
C2 = -c*(9/4); % Coeff. for O(gamma^2) term, found using Schwager-Poschel constant (valid for alp=bet=3/2)
end