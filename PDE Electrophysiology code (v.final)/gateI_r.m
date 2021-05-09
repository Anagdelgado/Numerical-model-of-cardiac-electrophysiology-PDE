function [gnou_r] = gateI_r(V_ant,g_r,dt)
ginf_r = (1+exp((20-V_ant)/6))^(-1);
tau_r = 9.5*exp(-(V_ant+40)^2/1800)+0.8;

gnou_r=(g_r+ginf_r*dt/tau_r)/(1+dt/tau_r);
%g = @(z)(z-g_r-(ginf_r-z)*dt/tau_r);
%gnou_r=fsolve(g,g_r,optimoptions('fsolve','Display','none'));