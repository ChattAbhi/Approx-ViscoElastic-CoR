function e=analyticCOR(alp,bet,gam,gtil,varargin)
% analyticCOR(alp,bet,gam,gtil) computes the $O(\gamma^2)$ Coefficient of
% Restitution, $e$ approximation. The input arguments alp, bet, gam, and 
% gtil, represent the rescaled parameters $\alpha$, $\beta$, $\gamma$, and 
% $\tilde{g}$, respectively. This computation approximates the integral 
% coefficients in the $O(\gamma^2)$ CoR expansion.
% Example:
% analyticCOR(3/2,3/2,1e-3,1e-3)
%
% analyticCOR(alp,bet,gam,gtil,choice,_) computes the $O(\gamma^2)$ 
% Coefficient of Restitution, $e$ approximation using the user-specified
% method indicated by the arguement choice. The argument choice can accept
% one of the two inputs 'approx' or 'numint'. If 'approx' is chosen the
% function computes the $O(\gamma^2)$ CoR expansion with the default
% approximated integral coefficient functions, whereas if 'numint' is 
% chosen the integral coefficient is numerically integrated. Additional
% settings such as absolute error tolerance, relative error tolerence, and 
% maximum number of iterations can be specified as Name-Value pair
% arguements, using the name keywords 'Abstol', 'RelTol', and 'MaxIter'
% respectively.
% Examples:
% analyticCOR(3/2,3/2,1e-3,1e-3,'approx')
% analyticCOR(3/2,3/2,1e-3,1e-3,'numint')
% analyticCOR(3/2,3/2,1e-3,1e-3,'nuint','AbsTol', 1e-6,'RelTol',1e-5,...
% 'MaxIter',500)

% Authors: Abhishek Chatterjee, Guillaume James, and Bernard Brogliato
% Address: Univ. Grenoble Alpes, INRIA, CNRS, Grenoble INP, LJK, Grenoble
%          38000 France 

choice='approx';
atol=1e-8; rtol=1e-7; m_iter=1000; 
if ~isempty(varargin)
    if ischar(varargin{1})
        if strcmpi(varargin{1},'approx') || strcmpi(varargin{1},'numint')
            choice=varargin{1};
            if strcmpi(varargin{1},'numint') && mod(length(varargin)-1,2)==0
                if ischar(varargin{2}) && ischar(varargin{4})
                    for i=2:2:length(varargin)-1
                        if strcmpi(varargin{i},'AbsTol')
                            if isnumeric(varargin{i+1})
                                atol=varargin{i+1};
                            else
                                error('Absolute Error Tolerence must be a numeric value!')
                            end
                        end
                        if strcmpi(varargin{i},'RelTol')
                            if isnumeric(varargin{i+1})
                                rtol=varargin{i+1};
                            else
                                error('Relative Error Tolerence must be a numeric value!')
                            end
                        end
                        if strcmpi(varargin{i},'MaxIter')
                            if isnumeric(varargin{i+1})
                                m_iter=varargin{i+1};
                            else
                                error('Maximum iteration must be a numeric value!')
                            end
                        end
                    end
                end
            end
        else
            error(['Unrecognized input ',varargin{1}])
        end
    else
        error('The choice must be a charcter array!')
    end
end
the=theta_gtilde(alp,gtil);
if strcmpi(choice,'approx')
    [I,Q]=VarCORCoeffsApprox(alp,bet,the);
elseif strcmpi(choice,'numint')
    [I,Q]=VarCORCoeffsNumeric(alp,bet,gtil,'AbsTol',atol,'RelTol',rtol,'MaxIter',m_iter);
end
%gam_c=2/(bet*(Q-(1/2)*I));

esqr=1 - 2*bet*I*gam + 2*(bet^2)*I*Q*(gam^2);
elin=1 - bet*I*gam + (bet^2)*I*(Q - (1/2)*I)*(gam^2);
e= min([sqrt(max([esqr,0])),max([elin,0])]);

% if gam>gam_c
%     esqr=1 - 2*bet*I*gam + 2*(bet^2)*I*Q*(gam^2);
%     e = sqrt(max([0,esqr])); 
% else
%     e = 1 - bet*I*gam + (bet^2)*I*(Q - (1/2)*I)*(gam^2);
% end

end