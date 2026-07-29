%% Project 1

close all
clear all
%% Inputs

Tx = 1;
Ty = 1;
l = 11;
m = 11;
delta_x = Tx/(m-1);%0.1;
delta_y = Tx/(m-1);%0.1;
alpha = 0.1;
Re = 30;
dP_dx = -12;
%%
yd = linspace(1,0,m)';
xd = linspace(0,1,l);
[u,v,C] = Main(Re,alpha,Tx,Ty,m,l);
%% Save Data

writematrix(u,'Data/uTable.csv')
writematrix(v,'Data/vTable.csv')
writematrix(C,'Data/CTable.csv')
%% Plots
% Developing Flow in a Channel

figure(Theme="light");
plot(yd,u)
legend("x* = "+string(xd))
ylabel("u")
xlabel("y*")
title("Developing Flow in a Channel")
saveas(gcf,"Figures/Developing_Flow.jpg")
quiver(xd,yd,u,v)
title("Developing Flow Vectors")
ylabel("y*")
xlabel("x*")
saveas(gcf,"Figures/Vectors.jpg")

ud_FD = u_FD(m);
u_err = zeros(size(u));
for i=1:l
    u_err(:,i) = (u(:,i)-ud_FD)./ud_FD;
end
u_err = rmmissing(u_err)*100;
figure(Theme="light");
hold on
plot(xd,max(u_err,[],1),'DisplayName',"Maximum")
plot(xd,mean(u_err,1),'DisplayName',"Average")
ylabel("% Error")
xlabel("x*")
legend
title("Developing Flow Error")
saveas(gcf,"Figures/FD_Error.jpg")
hold off

figure(Theme="light");
hold on
plot(xd(2:end),C)
ylabel("^{dp}/_{dx}")
title("Pressure Gradient")
xlabel("x*")
ylim([-25,0])
saveas(gcf,"Figures/Pressure_Gradient.jpg")
hold off
figure(Theme="light");
hold on
C_err = abs((C-dP_dx)/dP_dx)*100;
plot(xd(2:end),C_err)
ylabel("% Error")
title("Pressure Gradient Error")
xlabel("x*")
saveas(gcf,"Figures/Pressure_Gradient_Error.jpg")
hold off