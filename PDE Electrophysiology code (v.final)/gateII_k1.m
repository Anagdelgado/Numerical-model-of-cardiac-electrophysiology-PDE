function [gnou_k1] = gateII_k1(V_ant,phi_k)
constants

alpha_k1 = 0.1*(1+exp(0.06*(V_ant-phi_k-200)))^(-1);
beta_k1 = (3*exp(0.0002*(V_ant-phi_k+100))+exp(0.1*(V_ant-phi_k-10)))*(1+exp(-0.5*(V_ant-phi_k)))^(-1);

gnou_k1=alpha_k1/(alpha_k1+beta_k1);