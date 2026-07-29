%% Verification

clear all
load MainData.mat
%% Boundary Layer Polynomial Approximation

mid = (m+1)/2;
yd_m=yd(mid:end);
d_PA = 3/2*(1-1./u(mid,:));
eta = yd_m./d_PA(2:end);
eta(eta>1) = 1;
u_PA = 2*eta - eta.^2;

u_PA_full = [flip(u_PA,1);u_PA(2:end,:)];
figure(Theme="light");
plot(yd,u_PA_full)
title("Boundary Layer Polynomial Approximation")
ylabel("u")
xlabel("y*")
legend("x* = "+string(xd(2:end)))
%% Compare to Numerical Solution

u_err = rmmissing(abs(u(:,2:end)/u(mid,2:end)-u_PA_full)./u_PA_full)*100;
u_err_max = max(u_err,[],1);
u_err_mean = mean(u_err,1);

figure(Theme="light");
hold on
plot(xd(2:end),u_err_max,'DisplayName',"Maximum")
plot(xd(2:end),u_err_mean,'DisplayName',"Average")
title("Boundary Layer Polynomial Approximation Error")
ylabel("% Error")
xlabel("y*")
legend
hold off
%% Commencement Flow in Circular Pipe

dpdz=-12; %constant
v=0;
R=1;

% Roots of Bessel function of 1st kind
alpha = [2.4048, 5.5201, 8.6537, 11.7915, 14.9309];
% Source: https://mathworld.wolfram.com/BesselFunctionZeros.html

% Verify Roots
all(besselj(0,alpha)<1e-4)
u_BVP = @(t,r) sum( 8*besselj(0,alpha.*r/R) ./ (alpha.^3.*besselj(1,alpha)) );
u_BVP(0,R)<1e-4

r=linspace(0,1,mid);
u_CF = zeros(size(r))';
for ind=1:length(r)
    u_CF(ind) = u_BVP(0,r(ind));
end
u_CF = [flip(u_CF(2:end),1);u_CF];
figure(Theme="light");
plot(yd,u_CF)
title("Commencement Flow in Circular Pipe")
ylabel("u^{BVP}")
xlabel("y*")
u_err = zeros([m-2,l-1]);
for ind=2:l
u_err(:,ind-1) = (u(2:end-1,ind) - u_CF(2:end-1))./u_CF(2:end-1);
end
u_err = u_err*100;

%% Compare to Numerical Solution

u_err_max = max(u_err,[],1);
u_err_mean = mean(u_err,1);

figure(Theme="light");
hold on
plot(xd(2:end),u_err_max,'DisplayName',"Maximum")
plot(xd(2:end),u_err_mean,'DisplayName',"Average")
title("Commencement Flow Error")
ylabel("% Error")
xlabel("x*")
ylim([0,100])
legend
hold off