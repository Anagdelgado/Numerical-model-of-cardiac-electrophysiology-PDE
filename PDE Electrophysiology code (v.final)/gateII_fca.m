function [gnou_fca] = gateII_fca(V_ant,cant_ca,g_fca,dt)
ginf_fca = 0.685*((1+(cant_ca/0.000325)^8)^(-1)+0.1*(1+exp((cant_ca-0.0005)/0.0001))^(-1)+0.2*((1+exp((cant_ca-0.00075)/0.0008)))^(-1) +0.23);

if ginf_fca > g_fca && V_ant >= -60
    gnou_fca=g_fca;
else
    gnou_fca=(g_fca+ginf_fca*dt/2)/(1+dt/2);
    %g = @(z)(z-g_fca-(ginf_fca-z)*dt/2);
    %gnou_fca=fsolve(g,g_fca,optimoptions('fsolve','Display','off'));
end