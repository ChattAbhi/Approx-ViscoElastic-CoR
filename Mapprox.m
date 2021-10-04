function M=Mapprox(alp,bet,the,varargin)
% Mapprox(alp,bet,the) computes the analytical approximation of 
% $M_{\alpha,\beta}(\theta)$, with the inputs of alp, bet and the, which 
% represent the parameters $\alpha$, $\beta$, and $\theta$, respectively.

sel='Cheb';
if ~isempty(varargin)
    sel=varargin{1};
end

if strcmpi(sel,'pwise')
    if bet>=1 && bet<3/2
        M=Mfrac(alp,bet,the);
    elseif bet>=3/2 && bet<5/2
        M=Mxlnx(alp,bet,the);
    elseif bet>=5/2
        M=Mx2lnx(alp,bet,the);
    else
        error('Invalid!')
    end
elseif strcmpi(sel,'comb')
    kqr=10; krs=10; 
    M_q=Mfrac(alp,bet,the);
    if ~isfinite(M_q)
        M_q = 0;
    end
    M_r=Mxlnx(alp,bet,the);
    M_s=Mx2lnx(alp,bet,the);
    if ~isfinite(M_s)
        M_s = 0;
    end
    
    phi_qr = Phi(bet,kqr,1.4);
    phi_rs = Phi(bet,krs,2.4);
    
    Psi_q = 1 - phi_rs - phi_qr + phi_rs*phi_qr;
    Psi_r = phi_qr - phi_rs*phi_qr;
    Psi_s = phi_rs;
    
    
    M = Psi_q*M_q + Psi_r*M_r + Psi_s*M_s;
elseif strcmpi(sel,'all')
    M=Mall(alp,bet,the);
elseif strcmpi(sel,'frac')
    M=Mfrac(alp,bet,the);
elseif strcmpi(sel,'xlnx')
    M=Mxlnx(alp,bet,the);
elseif strcmpi(sel,'x2lnx')
    M=Mx2lnx(alp,bet,the);
elseif strcmpi(sel,'beta')
    M=((1-the)/(alp+the))*beta( (bet + the - 1/2)/(alp + the) + (the*(the-1))/(2*(alp + the)^2), 1/2 ) + (the/(alp+1))*beta(bet/(alp+1),1/2);
elseif strcmpi(sel,'cheb')
    %{
    npts=20;
    a=chebapprox(npts,alp,the);
    M=0;
    for k=0:npts
        M = M + (1/(alp+the))*a(k+1)*beta( (bet + k + the/2 - 1/2)/(alp+the)  ,1/2);
    end
    %}
    n=21;
    c=ChebCoeff(n,alp,the);
    A=zeros(n+1,n+1);
    v=(bet-0.5+the/2+(0:n))/(alp+the);
    A(1,:)=(1/(alp+the))*beta(v,1/2);
    A(2,1:n)=A(1,2:n+1) - A(1,1:n);
    for j=2:n
        A(j+1,1:n-j+1) = 2*A(j,2:n-j+2) - 2*A(j,1:n-j+1) - A(j-1,1:n-j+1);
    end
    %{
    B=zeros(n+1,n+1);
    for j=0:n
        for k=0:n-j
            B(j+1,k+1)=ajk(j,k,alp,bet,the);
        end
    end
    A=B;
    %}
    %{
    M=0;
    for j=0:n
        %M = M + c(j+1)*ajk(j,0,alp,bet,the);
        M = M + c(j+1)*A(j+1,1);
    end
    %}
    M=dot(c,A(:,1));
elseif strcmpi(sel,'ep')
    [c0,c1,c2,c3] = EpCoffs(alp,the);
    M = (1/(alp+the))*( c0*beta((bet + (c1*(1-the) + the)/2 - 1/2)/(alp+the),1/2) + c2*beta((bet + the/2 + 3/2)/(alp+the),1/2) + c3*beta( (bet + the/2 + 1/2)/(alp+the) , 1/2) );
else
    error('!!')
end
end

function M=Mfrac(alp,bet,the)
M0=(1/alp)*beta((2*bet-1)/(2*alp),1/2);
if bet>3/2
    dM0=(1/(2*alp^2))*( (alp-2*bet+1)*beta( (bet-1/2)/alp ,1/2 ) - (alp-2*bet+3)*beta( (bet-3/2)/alp,1/2 ) );
else
    dM0=nan;
end

M1=fun_M1(alp,bet,0);
dM1=fun_M1(alp,bet,1);
ddM1=fun_M1(alp,bet,2);

q0 = gamma(bet)*gamma(3/2-bet)/(sqrt(pi)*(1/2-bet));
q1 = M0 - M1 + dM1 + (3*q0)/2 - bet*q0;
q2 = 2*M1 - 2*M0 - dM1 - (5*q0)/2 + bet*q0;
q3 = M0;

M = q0*(the^(bet-1/2)) + q1*the^2 + q2*the + q3;
end

function M=Mxlnx(alp,bet,the)
M0=(1/alp)*beta((2*bet-1)/(2*alp),1/2);
if bet>3/2
    dM0=(1/(2*alp^2))*( (alp-2*bet+1)*beta( (bet-1/2)/alp ,1/2 ) - (alp-2*bet+3)*beta( (bet-3/2)/alp,1/2 ) );
else
    dM0=nan;
end

M1=fun_M1(alp,bet,0);
dM1=fun_M1(alp,bet,1);
ddM1=fun_M1(alp,bet,2);

r0 = 2*M0 - 2*M1 + 2*dM1 - ddM1;
r1 = M1 - M0 - dM1 + ddM1;
r2 = dM1 - ddM1;
r3 = M0;

M = r0*the*ln(the) + r1*the^2 + r2*the + r3;
end

function M=Mx2lnx(alp,bet,the)
M0=(1/alp)*beta((2*bet-1)/(2*alp),1/2);
if bet>3/2
    dM0=(1/(2*alp^2))*( (alp-2*bet+1)*beta( (bet-1/2)/alp ,1/2 ) - (alp-2*bet+3)*beta( (bet-3/2)/alp,1/2 ) );
else
    dM0=nan;
end

M1=fun_M1(alp,bet,0);
dM1=fun_M1(alp,bet,1);
ddM1=fun_M1(alp,bet,2);

s0 = 2*M0 - 2*M1 + dM0 + dM1;
s1 = M1 - M0 - dM0;
s2 = dM0;
s3=M0;

M=s0*(the^2)*ln(the) + s1*the^2 + s2*the + s3;
end

function M=Mall(alp,bet,the)
M0=(1/alp)*beta((2*bet-1)/(2*alp),1/2);
if bet>3/2
    dM0=(1/(2*alp^2))*( (alp-2*bet+1)*beta( (bet-1/2)/alp ,1/2 ) - (alp-2*bet+3)*beta( (bet-3/2)/alp,1/2 ) );
else
    dM0=nan;
end

M1=fun_M1(alp,bet,0);
dM1 = fun_M1(alp,bet,1);
ddM1 = fun_M1(alp,bet,2);
dddM1 = fun_M1(alp,bet,3);
ddddM1 = fun_M1(alp,bet,4);

b=[M0;M1;dM1;ddM1;dddM1;ddddM1];

A=[0,0,0,0,0,1; 1, 0, 0, 1 , 1 , 1; bet-1/2, 1, 1, 2, 1, 0;...
    (bet-1/2)*(bet-3/2), 1, 3, 2, 0, 0; ...
    (bet-1/2)*(bet-3/2)*(bet-5/2), -1, 2, 0, 0, 0;
    (bet-1/2)*(bet-3/2)*(bet-5/2)*(bet-7/2),2,-2, 0,0,0];
x=A\b;
M = x(1)*(the^(bet-1/2)) + x(2)*the*ln(the) + x(3)*(the^2)*ln(the) + x(4)*the^2 + x(5)*the + x(6);
end

function Mn0=fun_M0(alp,bet,n)
if n==0
    Mn0=(1/(alp+1))*beta(bet/(alp+1),1/2);
elseif n==1
    if bet>3/2
        Mn0=(1/(2*alp^2))*( (alp-2*bet+1)*beta( (bet-1/2)/alp ,1/2 ) - (alp-2*bet+3)*beta( (bet-3/2)/alp,1/2 ) );
    else
        Mn0=nan;
    end
else
    error('Cannot evaluate higher-order derivatives')
end
end

function Mn1=fun_M1(alp,bet,n)
P = ((((-1)^n)/(alp+1))*(prod(1:2:2*n-1)/(2^n)));
Q = 0;
for k=0:n
    cprod=1;
    for j=0:n-1
        cprod=cprod*( ( bet+k + (alp+1)*(1/2 - n + j) )/( (alp+1)*(1/2 - n + j) ) );
    end
    Q = Q + ((-1)^k)*nchoosek(n,k)*cprod*beta( (bet+k)/(alp+1), 1/2 );
end
Mn1=P*Q;
end

function phi=Phi(bet,k,s)
phi=1/(1+exp(-k*(bet-s)));
end

function a=chebapprox(npts,alp,the)
cheb_pt=[]; 
for k=0:npts
    xk = 0.5 + 0.5*cos( ((2*k+1)/(2*npts))*pi );
    yk = func(alp,the,xk);
    cheb_pt = [cheb_pt; [xk, yk]];
end
V=[]; Y=cheb_pt(:,2); 
for k=0:npts
    V=[V,cheb_pt(:,1).^k];
end

a=pinv(V)*Y;
end
function f=func(alp,the,x)
f = sqrt(((1-x^(alp+the))*(x^(1-the)))/(the + (1-the)*x - x^(alp+1)));
end

function [c0,c1,c2,c3] = EpCoffs(alp,the)
%EPCOFFS
%    [C0,C1,C2,C3] = EPCOFFS(ALP,THE)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    13-Jul-2021 03:58:08

t2 = alp+the;
t3 = alp.^2;
t4 = alp.^3;
t6 = alp.^5;
t8 = the.^2;
t9 = the.^3;
t11 = the.^5;
t13 = the.^7;
t14 = the-1.0;
t18 = the.^11;
t5 = t3.^2;
t7 = t3.^3;
t10 = t8.^2;
t12 = t8.^3;
t16 = t9.^3;
t17 = t8.^5;
t21 = 1.0./t2;
t22 = t4.*3.84e+2;
t23 = t11.*6.5e+1;
t24 = t11.*8.5e+1;
t26 = t9.*1.63e+2;
t28 = t9.*4.07e+2;
t31 = alp.*t11.*2.4e+1;
t32 = t4.*7.68e+2;
t36 = t4.*the.*2.4e+2;
t37 = t4.*the.*2.88e+2;
t38 = alp.*t9.*3.6e+2;
t41 = alp.*t8.*5.84e+2;
t42 = alp.*t9.*6.24e+2;
t43 = t3.*the.*7.6e+2;
t44 = alp.*t8.*1.42e+3;
t45 = t3.*the.*1.736e+3;
t49 = t3.*t11.*3.8e+1;
t53 = t3.*t8.*1.08e+2;
t54 = t3.*t8.*2.04e+2;
t55 = t4.*t8.*2.4e+2;
t56 = t3.*t9.*2.48e+2;
t57 = t4.*t8.*2.88e+2;
t58 = t3.*t9.*3.28e+2;
t15 = t10.^2;
t19 = t12.*9.0;
t20 = alp.*t12.*4.0;
t25 = t10.*1.47e+2;
t27 = t10.*2.67e+2;
t30 = t5.*the.*2.4e+1;
t34 = alp.*t10.*1.84e+2;
t35 = alp.*t10.*2.36e+2;
t39 = -t36;
t40 = -t37;
t46 = t5.*t8.*2.4e+1;
t47 = t3.*t10.*3.6e+1;
t48 = t5.*t9.*3.6e+1;
t51 = -t49;
t52 = t4.*t10.*9.6e+1;
t29 = -t20;
t33 = -t30;
t50 = -t48;
t59 = -t52;
t60 = t19+t22+t23+t25+t26+t31+t33+t34+t38+t39+t41+t43+t46+t47+t53+t55+t56;
t61 = t19+t24+t27+t28+t31+t32+t33+t35+t40+t42+t44+t45+t46+t47+t54+t57+t58;
t62 = 1.0./t60;
t63 = 1.0./t61;
c0 = (t62.*t63.*(t7.*8.84736e+5+t12.*2.26981e+5+t13.*3.3489e+5+t15.*2.20515e+5+t16.*8.19e+4+t17.*1.8075e+4+t18.*2.25e+3+alp.*t11.*1.652124e+6+alp.*t12.*2.026908e+6+alp.*t13.*1.155024e+6+alp.*t15.*3.8412e+5+alp.*t16.*7.974e+4+alp.*t17.*9.9e+3+alp.*t18.*6.0e+2+t3.*t10.*5.08008e+6+t4.*t9.*8.44192e+6+t5.*t8.*7.99488e+6+t3.*t11.*4.841532e+6+t4.*t10.*5.537952e+6+t5.*t9.*2.688192e+6-t6.*t8.*2.7648e+4+t3.*t12.*2.395188e+6+t4.*t11.*2.602752e+6+t5.*t10.*1.839024e+6+t6.*t9.*1.059264e+6+t7.*t8.*3.73248e+5+t3.*t13.*7.58808e+5+t4.*t12.*8.55936e+5+t5.*t11.*5.50656e+5+t6.*t10.*8.1216e+4-t7.*t9.*8.4672e+4+t4.*t13.*2.18496e+5+t5.*t12.*1.66752e+5+t6.*t11.*9.1584e+4+t7.*t10.*4.6656e+4+t3.*t15.*1.69272e+5+t5.*t13.*2.7072e+4+t6.*t12.*8.64e+3-t7.*t11.*5.184e+3+t3.*t16.*2.43e+4+t4.*t15.*3.4272e+4+t6.*t13.*3.456e+3+t7.*t12.*1.728e+3+t3.*t17.*1.86e+3+t4.*t16.*3.392e+3+t5.*t15.*4.464e+3+t6.*the.*4.091904e+6-t7.*the.*3.31776e+5+t10.^3.*1.25e+2))./3.0;
if nargout > 1
    c1 = (t21.*(t9.*2.7e+1+t10.*1.8e+1+t11.*3.0+alp.*t8.*8.4e+1+alp.*t9.*5.2e+1+alp.*t10.*8.0+t3.*t8.*6.8e+1+t3.*t9.*1.2e+1+t3.*the.*7.2e+1+t4.*the.*6.4e+1+t5.*the.*8.0).*(3.0./2.0))./(t3.*9.6e+1+t8.*6.1e+1+t9.*3.0e+1+t10.*5.0+alp.*t8.*3.6e+1+alp.*t9.*8.0+alp.*the.*1.48e+2+t3.*t8.*1.2e+1-t3.*the.*1.2e+1);
end
if nargout > 2
    c2 = t14.*t21.*t63.*(t10.*2.9e+1+t11.*4.7e+1+t12.*1.9e+1+t13+t29+t50+t51+t59+alp.*t9.*1.96e+2+alp.*t10.*2.92e+2+alp.*t11.*9.2e+1+t3.*t8.*4.76e+2+t3.*t9.*6.86e+2+t4.*t8.*9.36e+2+t3.*t10.*1.72e+2+t4.*t9.*2.16e+2+t5.*t8.*1.8e+2+t4.*the.*3.36e+2+t5.*the.*5.76e+2+t6.*the.*1.44e+2).*(-1.0./6.0);
end
if nargout > 3
    c3 = (t14.*t21.*t62.*(t10.*1.49e+2+t11.*1.91e+2+t12.*4.3e+1+t13+t29+t50+t51+t59+alp.*t9.*9.4e+2+alp.*t10.*1.156e+3+alp.*t11.*2.12e+2+t3.*t8.*1.82e+3+t3.*t9.*2.702e+3+t4.*t8.*3.384e+3+t3.*t10.*4.12e+2+t4.*t9.*3.6e+2+t5.*t8.*3.24e+2+t4.*the.*1.056e+3+t5.*the.*1.728e+3+t6.*the.*2.88e+2))./6.0;
end
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

function a=ajk(j,k,alp,bet,the)
if j==0
    a=(1/(alp+the))*beta((bet + k - 1/2 + the/2)/(alp + the),1/2);
elseif j==1
    a=ajk(j-1,k+1,alp,bet,the) - ajk(j-1,k,alp,bet,the);
else
    a=2*ajk(j-1,k+1,alp,bet,the) - 2*ajk(j-1,k,alp,bet,the) - ajk(j-2,k,alp,bet,the);
end
end
