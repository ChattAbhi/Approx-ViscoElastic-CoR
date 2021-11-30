function e=getCOR(m,k,gam0,v0,alp,bet,g,F,varargin)
% getCOR(m,k,gam0,v0,alp,bet,g,F,method) computes the CoR using the physical
% parameters of m (mass), k (stiffness constant), gam0 (damping ratio, ie 
% c/k, where c is the damping coefficient), v0 (positive pre-impact 
% velocity, direction is assumed positive towards the contact >0), alp 
% ($\alpha$ exponent on the stiffness term >=1), bet ($\beta$ exponent on 
% the damping term >=1), g (gravitational acceleration constant), and F 
% (constant compressive force value >0). The argument method can be used
% to select the method for computation of the CoR. The list of allowed
% method inputs and their respective descriptions are provided below.
%
% 
% * 'ord2approx'   -  (Default) Computes CoR using the $O(\gamma^2)$ 
%                     expansion, while using approximate analytical
%                     solution for the integral coefficient functions.
%
% * 'ord1approx'   -  Computes the CoR using the $O(\gamma)$ (first-order)
%                     approximation 
% 
% * 'ord2schw-pos' -  Computes the CoR when $\beta = 3/2$ with a 
%                     $O(\gamma^2)$ expansion using the second-order 
%                     coefficient reported by Schwager and Poschel in
%                     "Coefficient of restitution of colliding viscoelastic
%                     spheres: The effect of delayed recovery", Physical
%                     Review E (2008)
%
% * 'ord2numint'   -  Computes the CoR using the $O(\gamma^2)$ expansion,
%                     with coefficient function computed by means of
%                     numerical integration
%
% * 'dirnumint'    -  Computes the CoR by numerically integration the
%                     dynamic equations given by the viscoelastic contact 
%                     model
% Examples:
% getCOR(1.54e-1,3.6138e10,1.5237e-6,0.1,3/2,3/2,9.8,0)
% getCOR(1.54e-1,3.6138e10,1.5237e-6,0.1,3/2,3/2,9.8,0,'ord2approx')
% getCOR(1.54e-1,3.6138e10,1.5237e-6,0.1,3/2,3/2,9.8,0,'ord1approx')
% getCOR(1.54e-1,3.6138e10,1.5237e-6,0.1,3/2,3/2,9.8,0,'ord2schw-pos')
% getCOR(1.54e-1,3.6138e10,1.5237e-6,0.1,3/2,3/2,9.8,0,'ord2numint')
% getCOR(1.54e-1,3.6138e10,1.5237e-6,0.1,3/2,3/2,9.8,0,'dirnumint')
%
%
% getCOR(m,k,gam0,v0,alp,bet,g,F,method,_) computes the CoR for a
% viscoelastic impact with additional settings passed in as Name-Value
% pair arguments, for the methods using numerical integation namely
% 'ord2numint' and 'dirnumint'. Valid name keywords include 'AbsTol',
% 'RelTol', and 'MaxIter', which can be used to set the absolute error
% tolerence, relative error tolerence, and the maximum number of iterations
% for the numerical integration. Additionally when using the method 
% 'dirnumint', a tolerence for the convergence to the equilibrium point
% can be set using the name keyword 'EquTol'. 
% Examples:
% getCOR(1.54e-1,3.6138e10,1.5237e-6,0.1,3/2,3/2,9.8,0,'ord2numint',...
% 'AbsTol',1e-6,'RelTol',1e-5,'MaxIter',100)
% getCOR(1.54e-1,3.6138e10,1.5237e-6,0.1,3/2,3/2,9.8,0,'dirnumint',...
% 'AbsTol',1e-6,'RelTol',1e-5,'MaxIter',100,'EquTol',1e-5)

% Authors: Abhishek Chatterjee, Guillaume James, and Bernard Brogliato
% Address: Univ. Grenoble Alpes, INRIA, CNRS, Grenoble INP, LJK, Grenoble
%          38000 France 

choice='ord2approx';

if ~isempty(varargin)
    if ischar(varargin{1})
        if ~isempty(find(strcmpi({'ord2approx','ord2numint','ord2schw-pos','ord1approx','dirnumint'},varargin{1}),1))
            choice=varargin{1};
        else
            error(['The specified method ',varargin{1},' is invalid!'])
        end
    else
        error('The method specified must be a character array!');
    end
    if length(varargin)>1
        if  mod(length(varargin)-1,2)==0
            for i=2:2:length(varargin)-1
                if ischar(varargin{i})
                    if isempty(find(strcmpi({'AbsTol','RelTol','MaxIter','EquTol'},varargin{i}),1))
                        error(['The name keyword ',varargin{i},' is unrecognized!']);
                    end
                else
                    error('Name keyword must be a character array!')
                end
            end
        else
            error('Uneven Name-Value pairs!')
        end
    end
end

[gam,gtil]=ScaledParams(m,k,gam0,v0,alp,bet,g,F);

if strcmpi(choice,'ord2approx')
    e=analyticCOR(alp,bet,gam,gtil);
elseif strcmpi(choice,'ord2numint')
    atol=1e-8; rtol=1e-7; m_iter=100; 
    if length(varargin)>1
        for i=2:2:length(varargin)-1
            if strcmpi(varargin{i},'AbsTol')
                atol=varargin{i+1};
            end
            if strcmpi(varargin{i},'RelTol')
                rtol=varargin{i+1};
            end
            if strcmpi(varargin{i},'MaxIter')
                m_iter=varargin{i+1};
            end
        end
    end
    e=analyticCOR(alp,bet,gam,gtil,'numint','AbsTol',atol,'RelTol',rtol,'MaxIter',m_iter);
elseif strcmpi(choice,'ord2schw-pos')
    if bet==3/2
        [C0,C1,C2]=ConstCORCoeffs(alp,bet);
        e = max([1 - gam*C0 - gam*gtil*C1 - C2*gam^2,0]);
    else
        error('Schwager-Poschel method is only valid when $\beta=3/2$!')
    end
elseif strcmpi(choice,'ord1approx')
    [C0,C1]=ConstCORCoeffs(alp,bet);
    e = max([1 - gam*C0 - gam*gtil*C1,0]);
elseif strcmpi(choice,'dirnumint')
    atol=1e-14; rtol=1e-10; etol=1e-8; m_iter=100;
    if length(varargin)>1
        for i=2:2:length(varargin)-1
            if strcmpi(varargin{i},'AbsTol')
                atol=varargin{i+1};
            end
            if strcmpi(varargin{i},'RelTol')
                rtol=varargin{i+1};
            end 
            if strcmpi(varargin{i},'MaxIter')
                m_iter=varargin{i+1};
            end
            if strcmpi(varargin{i},'EquTol')
                etol=varargin{i+1};
            end
        end
    end
    e=numericCOR(alp,bet,gam,gtil,'AbsTol',atol,'RelTol',rtol,'MaxIter',m_iter,'Equtol',etol);
else
    error(['The specified method ',choice,' is invalid!'])
end

end