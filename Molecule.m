function [u,v,coeff] = Molecule(dx,dy,Re,alpha,uw,vw,up,vp)
%% MOLECULE Computational Molecule
arguments (Input)
    dx
    dy
    Re
    alpha
    uw
    vw
    up
    vp
end
arguments (Output)
    u % x velocity component
    v % y velocity component
    coeff
end
m = 1/dy+1;
%% 
% $$A_p u_p =A_S u_S +A_N u_N +B_u +C$$
%% Coefficients
% $$A_N =-\frac{1}{\Delta y^2 }+\frac{{\mathrm{Re}}_H \alpha v}{\Delta y}$$
An = -1/dy^2 + Re*alpha*vp/dy;
%% 
% $$A_S =-\frac{1}{\Delta y^2 }-\frac{{\mathrm{Re}}_H \alpha v}{\Delta y}$$
As = -1/dy^2 - Re*alpha*vp/dy;
Ap = An + As;
%% 
% $$B_u ={\mathrm{Re}}_H \;\alpha \;u_p \frac{\;u_p -u_w }{\Delta x}-\alpha^2 
% \frac{\;u_w +u_p -2u_p }{\Delta x^2 }$$
Bu = Re * alpha * up .* (up - uw)/dx - alpha^2 .* (uw + up - 2*up)/dx^2;
%% 
% Iterate to find $\frac{\mathit{dp}}{\mathit{dx}}$
mass = @(x) sum(TDMA(An,As,Ap,Bu,x,m)*dy)-1;
C = fzero(mass,0);
%% 
% Velocity Components
[u,P,Q] = TDMA(An,As,Ap,Bu,C,m);
v = (u-up)/dx*dy .* [0;ones([(m-3)/2 1]);0;-ones([(m-3)/2 1]);0];
mid = (m-1)/2;
v(mid)=mean([v(mid-1),v(mid+1)]);
%% 
% Save Coefficients
coeff.C  = C;
coeff.An = An;
coeff.As = As;
coeff.Ap = Ap;
coeff.Bu = Bu;
coeff.P = P;
coeff.Q = Q;
end
%% TDMA Routine
function [u,P,Q] = TDMA(An,As,Ap,Bu,C,m)
arguments (Input)
    An
    As
    Ap
    Bu
    C
    m
end
arguments (Output)
    u
    P
    Q
end
D = Bu + C;
u = zeros([m,1]);
P = zeros([m,1]);
Q = zeros([m,1]);
for j=m-1:-1:2
    a = Ap(j) - As(j)*P(j+1);
    b = An(j);
    d = D(j) + As(j)*Q(j+1);
    P(j) = b/a;
    Q(j) = d/a;
end
% Update velocity field for the next iteration
for j=linspace(2, m-1, m-2)
    u(j) = P(j) * u(j-1) + Q(j);
end
end