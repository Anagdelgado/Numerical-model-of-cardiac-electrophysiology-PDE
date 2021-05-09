function [gnou_xs] = gateI_xs(V_ant,g_xs,dt)
ginf_xs = (1+exp((-5-V_ant)/14))^(-1);
tau_xs = 1100*(1+exp((-10-V_ant)/6))^(-1/2)*(1+exp((V_ant-60)/20))^(-1);

gnou_xs=(g_xs+ginf_xs*dt/tau_xs)/(1+dt/tau_xs);
%g = @(z)(z-g_xs-(ginf_xs-z)*dt/tau_xs);
%gnou_xs=fsolve(g,g_xs,optimoptions('fsolve','Display','off'));