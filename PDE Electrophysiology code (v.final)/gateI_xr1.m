function [gnou_xr1] = gateI_xr1(V_ant,g_xr1,dt)
ginf_xr1 = (1+exp((-26-V_ant)/7))^(-1);
tau_xr1 = 2700*(1+exp((-45-V_ant)/10))^(-1)*((1+exp((V_ant+30)/11.5)))^(-1);

gnou_xr1=(g_xr1+ginf_xr1*dt/tau_xr1)/(1+dt/tau_xr1);
%g = @(z)(z-g_xr1-(ginf_xr1-z)*dt/tau_xr1);
%gnou_xr1=fsolve(g,g_xr1,optimoptions('fsolve','Display','off'));