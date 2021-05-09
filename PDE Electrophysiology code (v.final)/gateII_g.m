function [gnou_g] = gateII_g(V_ant,cant_ca,g_g,dt)
if cant_ca <= 0.00035
    ginf_g = (1+(cant_ca/0.00035)^6)^(-1);
else
    ginf_g = (1+(cant_ca/0.00035)^16)^(-1);
end

if ginf_g > g_g && V_ant >= -60
    gnou_g=g_g;
else
    gnou_g=(g_g+ginf_g*dt/2)/(1+dt/2);
    %g = @(z)(z-g_g-(ginf_g-z)*dt/2);
    %gnou_g=fsolve(g,g_g,optimoptions('fsolve','Display','off'));
end