function um=maxcompress_zero(alp)
% maxcompress_zero(alp) computes the value of maximum compression
% corresponding to the case $\tilde{g}=0$. The input argument alp represent
% the parameter $\alpha$.

% Authors: Abhishek Chatterjee, Guillaume James, and Bernard Brogliato
% Address: Univ. Grenoble Alpes, INRIA, CNRS, Grenoble INP, LJK, Grenoble
%          38000 France 
um=((alp+1)/2)^(1/(alp+1));
end