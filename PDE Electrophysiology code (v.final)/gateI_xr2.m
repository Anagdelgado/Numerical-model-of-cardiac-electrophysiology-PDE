function [gnou_xr2] = gateI_xr2(V_ant,g_xr2,dt)
ginf_xr2 = (1+exp((V_ant+88)/24))^(-1);
tau_xr2 = 3.36*(1+exp((-60-V_ant)/20))^(-1)*(1+exp((V_ant-60)/20))^(-1);

gnou_xr2=(g_xr2+ginf_xr2*dt/tau_xr2)/(1+dt/tau_xr2);
%g = @(z)(z-g_xr2-(ginf_xr2-z)*dt/tau_xr2);
%gnou_xr2=fsolve(g,g_xr2,optimoptions('fsolve','Display','off'));