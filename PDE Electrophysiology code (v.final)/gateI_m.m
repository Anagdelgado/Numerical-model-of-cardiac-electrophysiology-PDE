function [gnou_m] = gateI_m(V_ant,g_m,dt)
ginf_m = (1+exp((-56.86-V_ant)/9.03))^(-2);
tau_m = 0.1*(1+exp((-60-V_ant)/5))^(-1)*((1+exp((V_ant+35)/5))^(-1)+(1+exp((V_ant-50)/200))^(-1));


gnou_m=(g_m+ginf_m*dt/tau_m)/(1+dt/tau_m);
%g = @(z)(z-g_m-(ginf_m-z)*dt/tau_m);
%gnou_m=fsolve(g,g_m,optimoptions('fsolve','Display','off'));
