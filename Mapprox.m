function M=Mapprox(alp,bet,the,varargin)
% Mapprox(alp,bet,the) computes the analytical approximation of 
% $M_{\alpha,\beta}(\theta)$, with the inputs of alp, bet and the, which 
% represent the parameters $\alpha$, $\beta$, and $\theta$, respectively.

% Authors: Abhishek Chatterjee, Guillaume James, and Bernard Brogliato
% Address: Univ. Grenoble Alpes, INRIA, CNRS, Grenoble INP, LJK, Grenoble
%          38000 France 

n=21;
c=ChebCoeff(n,alp,the);
A=zeros(n+1,n+1);
v=(bet-0.5+the/2+(0:n))/(alp+the);
A(1,:)=(1/(alp+the))*beta(v,1/2);
A(2,1:n)=A(1,2:n+1) - A(1,1:n);
for j=2:n
    A(j+1,1:n-j+1) = 2*A(j,2:n-j+2) - 2*A(j,1:n-j+1) - A(j-1,1:n-j+1);
end
M=dot(c,A(:,1));

end

function f=func(alp,the,x)
f = sqrt(((1-x^(alp+the))*(x^(1-the)))/(the + (1-the)*x - x^(alp+1)));
end

function c=ChebCoeff(n,alp,the)
c = zeros(1,n+1);
for j=0:n
    c(j+1) = 0;
    for k = (n+1)/2:n
        xk = cos(pi*(2*k+1)/(2*n+2)) + 1;
        if j==0
            c(j+1) = c(j+1) + (2/(n+1))*func(alp,the,xk);
        elseif j>=2 && mod(j,2)==0
            c(j+1) = c(j+1) + (4/(n+1))*func(alp,the,xk)*cos(j*pi*(2*k+1)/(2*n+2));
        end
    end
end
end