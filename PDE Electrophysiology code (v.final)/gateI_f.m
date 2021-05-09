function [gnou_f] = gateI_f(V_ant,g_f,dt)
ginf_f = (1+exp((V_ant+20)/7))^(-1);
tau_f = 1125*exp(-(V_ant+27)^2/240)+165*(1+exp((25-V_ant)/10))^(-1)+80;

gnou_f=(g_f+ginf_f*dt/tau_f)/(1+dt/tau_f);
%g = @(z)(z-g_f-(ginf_f-z)*dt/tau_f);
%gnou_f=fsolve(g,g_f,optimoptions('fsolve','Display','none'));