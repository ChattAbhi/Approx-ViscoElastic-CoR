function [gam,gtil]=ScaledParams(m,k,gam0,v0,alp,bet,g,F)
% [gam,gtil]=ScaledParams(m,k,gam0,alp,bet,g,F) converts the physical
% dimensional parameters of m (mass), k (stiffness), gam0 
% (damping ratio: c/k, where c is the damping coefficient), v0 (positive 
% pre-impact velocity), alp ($\alpha$ is the exponent on the stiffness 
% term), bet ($\beta$ is the exponent on the damping term), 
% g (gravitational acceleration constant), and F (compressive external 
% force), using appropriate units, to the scaled non-dimensional parameters 
% of gam ($\gamma$) and gtil ($\tilde{g}$).
% Example:
% [gam,gtil]=ScaledParams(1.54e-1,3.6138e10,1.5237e-6,0.1,3/2,3/2,9.8,0)

% Authors: Abhishek Chatterjee, Guillaume James, and Bernard Brogliato
% Address: Univ. Grenoble Alpes, INRIA, CNRS, Grenoble INP, LJK, Grenoble
%          38000 France 

if numel(m)~=1 || numel(k)~=1 || numel(v0)~=1 || numel(gam0)~=1 || numel(alp)~=1 || numel(bet)~=1 || numel(g)~=1 || numel(F)~=1
    error('All physical parameter inputs must be scalar non-negative numeric values!')
end
if ~isnumeric([m,k,gam0,v0,alp,bet,g,F])
    error('All physical parameter inputs must be non-negative numeric values!')
end
if ~isempty(find([m,k,gam0,v0,alp,bet,g,F]<0,1))
    error('All physical parameter inputs must be non-negative numeric values!')
end
if v0<=0 || m<=0
    error('The pre-impact velocity and mass mus be strictly positive values!')
end
if alp<1 || bet<1
    error('The parameters alp and bet must have values greater or equal to 1!')
end

gam=gam0*(v0^((2*bet)/(alp+1) - 1))*((k/m)^(1-(bet/(alp+1))));
gtil=((m/k)^(1/(alp+1)))*(v0^(-2*alp/(alp+1)))*(g+F/m);

end