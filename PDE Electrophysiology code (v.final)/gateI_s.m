function [gnou_s] = gateI_s(V_ant,g_s,dt)

% Epicardium
ginf_s = (1+exp((V_ant+20)/5))^(-1);
tau_s = 85*exp(-(V_ant+45)^2/320)+5*(1+exp((V_ant-20)/5))^(-1)+3;

% Endocardium
%ginf_s = 1+exp((V_ant+28)/5);
%tau_s = 1000*exp(-(V_ant+67)^2/1000)+8;

gnou_s=(g_s+ginf_s*dt/tau_s)/(1+dt/tau_s);
%g = @(z)(z-g_s-(ginf_s-z)*dt/tau_s);
%gnou_s=fsolve(g,g_s,optimoptions('fsolve','Display','none'));