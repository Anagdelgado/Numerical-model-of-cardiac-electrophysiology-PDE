function [gnou_h] = gateI_h(V_ant,g_h,dt)
ginf_h = (1+exp((V_ant+71.55)/7.43))^(-2);
if (V_ant >= -40)
    tau_h = 0.1688*(1+exp(-(V_ant+10.66)/11.1));
elseif (V_ant < -40)
    tau_h = (0.057*exp(-(80+V_ant)/6.8)+2.7*exp(0.079*V_ant)+3.1*(10^5)*exp(0.3485*V_ant))^(-1);
end

gnou_h=(g_h+ginf_h*dt/tau_h)/(1+dt/tau_h);
%g = @(z)(z-g_h-(ginf_h-z)*dt/tau_h);
%gnou_h=fsolve(g,g_h,optimoptions('fsolve','Display','none'));