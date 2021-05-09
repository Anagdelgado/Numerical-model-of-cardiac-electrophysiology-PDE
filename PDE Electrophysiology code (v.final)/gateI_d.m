function [gnou_d] = gateI_d(V_ant,g_d,dt)
ginf_d = (1+exp((-5-V_ant)/7.5))^(-1);
tau_d = (1.4*(1+exp((-35-V_ant)/13))^(-1)+0.25)*(1.4*(1+exp((V_ant+5)/5))^(-1))+(1+exp((50-V_ant)/20))^(-1);

gnou_d=(g_d+ginf_d*dt/tau_d)/(1+dt/tau_d);
%g = @(z)(z-g_d-(ginf_d-z)*dt/tau_d);
%gnou_d=fsolve(g,g_d,optimoptions('fsolve','Display','off'));