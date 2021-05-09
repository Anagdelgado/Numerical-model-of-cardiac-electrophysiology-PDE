function [phi_ionx] = phi_ion(R,T,F,z_ion,c_ion0,c_ion)
phi_ionx = R*T*log(c_ion0/c_ion)/(z_ion*F);