function [gnou_j] = gateI_j(V_ant,g_j,dt)
ginf_j = (1+exp((V_ant+71.55)/7.43))^(-2);
if (V_ant >= -40)
    a_j = 0;
    b_j = 0.6*exp(0.057*V_ant)*(1+exp(-0.1*(V_ant+32)))^(-1);
end
if (V_ant < -40)
    a_j = (-2.5428*(10^4)*exp(0.2444*V_ant)-6.948*(10^(-6))*exp(-0.04391*V_ant))*(V_ant+37.78)*(1+exp(0.311*(V_ant+79.23)))^(-1);
    b_j = 0.02424*exp(-0.01052*V_ant)*(1+exp(-0.1378*(V_ant+40.14)))^(-1);
end
tau_j = (a_j+b_j)^(-1);

gnou_j=(g_j+ginf_j*dt/tau_j)/(1+dt/tau_j);
%g = @(z)(z-g_j-(ginf_j-z)*dt/tau_j);
%gnou_j=fsolve(g,g_j,optimoptions('fsolve','Display','none'));