function [vect_x,vect_r,iter_fi,gnou_k1,gnou_fca,gnou_g] = newton_iter_local(niter,tol,cant_k,cant_na,cant_ca,cant_srca,...
            gnou_m,gnou_h,gnou_j,gnou_xr1,gnou_xr2,gnou_xs,gnou_r,gnou_s,gnou_d,gnou_f,gant_fca,gant_g,V_ant,dt)

constants

% Functions for phi_ion, phi_ks, y_ca, y_srca
phi_k = @(z1) phi_ion(RC,T,F,z_k,c_k0,z1);
phi_na = @(z2) phi_ion(RC,T,F,z_na,c_na0,z2);
phi_ca = @(z3) phi_ion(RC,T,F,z_ca,c_ca0,z3);
phi_ks = @(z1,z2) RC*T/F*log((c_k0+p_kna*c_na0)*(z1+p_kna*z2)^(-1));
y_ca = @(z3) (1+(c_tot*c_buf)*(z3+c_buf)^(-2))^(-1);
y_srca = @(z4) (1+(c_srtot*c_srbuf)*(z4+c_srbuf)^(-2))^(-1);

% Functions for gateII gate variables (3)
f_gnou_k1 = @(z1) gateII_k1(V_ant,phi_k(z1));
f_gnou_fca = @(z3) gateII_fca(V_ant,z3,gant_fca,dt);
f_gnou_g = @(z3) gateII_g(V_ant,z3,gant_g,dt);

% Functions for currents
I_na = @(z2) Cmax_na*gnou_m^3*gnou_h*gnou_j*(V_ant-phi_na(z2));
I_bna = @(z2) Cmax_bna*(V_ant-phi_na(z2));
I_nak = @(z2) Imax_nak*(c_k0*z2)*((z2+c_nak)*(c_k0+c_kna)*(1+0.1245*exp(-0.1*V_ant*F/(RC*T))+0.0353*exp(-V_ant*F/(RC*T))))^(-1);
I_naca = @(z2,z3) Imax_naca*(exp(y*V_ant*F/(RC*T))*z2^3*c_ca0-exp((y-1)*V_ant*F/(RC*T))*c_na0^3*z3*y_naca)...
            *((c_naca^3+c_na0^3)*(c_cana+c_ca0)*(1+k_naca*exp((y-1)*V_ant*F/(RC*T))))^(-1);
I_k1 = @(z1) Cmax_k1*f_gnou_k1(z1)*(c_k0/5.4)^(1/2)*(V_ant-phi_k(z1));
I_kr = @(z1) Cmax_kr*gnou_xr1*gnou_xr2*(c_k0/5.4)^(1/2)*(V_ant-phi_k(z1));
I_ks = @(z1,z2) Cmax_ksepi*gnou_xs^2*(V_ant-phi_ks(z1,z2));
I_pk = @(z1) Cmax_pk*(1+exp((25-V_ant)/5.98))^(-1)*(V_ant-phi_k(z1));
I_t0 = @(z1) Cmax_t0epi*gnou_r*gnou_s*(V_ant-phi_k(z1));
I_cal = @(z3) Cmax_cal*gnou_d*gnou_f*f_gnou_fca(z3)*4*F^2*V_ant*(RC*T)^(-1)*(z3*exp(2*V_ant*F*(RC*T)^(-1))-0.341*c_ca0)*(exp(2*V_ant*F*(RC*T)^(-1))-1)^(-1);
I_bca = @(z3) Cmax_bca*(V_ant-phi_ca(z3));
I_pca = @(z3) Cmax_pca*z3*(c_pca+z3)^(-1);
I_leak = @(z3,z4) Imax_leak*(z4-z3);
I_up = @(z3) Imax_up*(1+c_up^2/z3^2)^(-1);
I_rel = @(z3,z4) Imax_rel*gnou_d*f_gnou_g(z3)*(1+y_rel*z4^2*(c_rel^2+z4^2)^(-1));
        
% Functions for concentration residuals (implicit)
f1_k = @(z1,z2) (z1 - cant_k + C*(Vol*F)^(-1)*(I_k1(z1)+I_kr(z1)+I_ks(z1,z2)-2*I_nak(z2)+I_pk(z1)+I_t0(z1))*dt);
f2_na = @(z2,z3) (z2 - cant_na + C*(Vol*F)^(-1)*(I_na(z2)+I_bna(z2)+3*I_nak(z2)+3*I_naca(z2,z3))*dt);
f3_ca = @(z2,z3,z4) (z3 - cant_ca + (C*(2*Vol*F)^(-1)*(I_cal(z3)+I_bca(z3)+I_pca(z3)-2*I_naca(z2,z3))-I_leak(z3,z4)+I_up(z3)-I_rel(z3,z4))*y_ca(z3)*dt);
f4_srca = @(z3,z4) (z4 - cant_srca + Vol*Vol_sr^(-1)*(I_leak(z3,z4)-I_up(z3)+I_rel(z3,z4))*y_srca(z4)*dt);


% 4 concentrations residuals vector
RES = @(z1,z2,z3,z4)[f1_k(z1,z2) f2_na(z2,z3) f3_ca(z2,z3,z4) f4_srca(z3,z4)];

% Derivatives of current functions
dk_I_k1 = @(z1) Cmax_k1*f_gnou_k1(z1)*(c_k0/5.4)^(1/2)*RC*T*(z_k*F*z1)^(-1);
dk_I_kr = @(z1) Cmax_kr*gnou_xr1*gnou_xr2*(c_k0/5.4)^(1/2)*RC*T*(z_k*F*z1)^(-1);
dk_I_ks = @(z1,z2) Cmax_ksepi*gnou_xs^2*RC*T*(F*(z1+ p_kna*z2))^(-1);
dk_I_pk = @(z1) Cmax_pk*(1+exp((25-V_ant)/5.98))^(-1)*RC*T*(z_k*F*z1)^(-1);
dk_I_t0 = @(z1) Cmax_t0epi*gnou_r*gnou_s*RC*T*(z_k*F*z1)^(-1);

dna_I_ks = @(z1,z2) Cmax_ksepi*gnou_xs^2*RC*T*p_kna*(F*(z1+ p_kna*z2))^(-1);
dna_I_na = @(z2) Cmax_na*gnou_m^3*gnou_h*gnou_j*RC*T*(z_na*F*z2)^(-1);
dna_I_bna = @(z2) Cmax_bna*RC*T*(z_na*F*z2)^(-1);
dna_I_nak = @(z2) Imax_nak*c_k0*c_nak*((c_k0+c_kna)*(1+0.1245*exp(-0.1*V_ant*F/(RC*T))+0.0353*exp(-V_ant*F/(RC*T))))^(-1)*(z2+c_nak)^(-2);
dna_I_naca = @(z2) 3*Imax_naca*exp(y*V_ant*F/(RC*T))*z2^2*c_ca0*((c_naca^3+c_na0^3)*(c_cana+c_ca0)*(1+k_naca*exp((y-1)*V_ant*F/(RC*T))))^(-1);

dca_I_cal = @(z3) Cmax_cal*gnou_d*gnou_f*f_gnou_fca(z3)*4*F^2*V_ant*(RC*T)^(-1)*exp(2*V_ant*F*(RC*T)^(-1))*(exp(2*V_ant*F*(RC*T)^(-1)))^(-1);
dca_I_naca = (-1)*Imax_naca*(exp((y-1)*V_ant*F/(RC*T))*c_na0^3*y_naca)*((c_naca^3+c_na0^3)*(c_cana+c_ca0)*(1+k_naca*exp((y-1)*V_ant*F/(RC*T))))^(-1);
dca_I_bca = @(z3) Cmax_bca*RC*T*(z_ca*F*z3)^(-1);
dca_I_pca = @(z3) Cmax_pca*c_pca*(c_pca+z3)^(-2);
dca_I_leak = (-1)*Imax_leak;
dca_I_up = @(z3) 2*Imax_up*c_up^2*z3*(c_up^2+z3^2)^(-2);

dsrca_I_leak = Imax_leak;
dsrca_I_rel = @(z3,z4) 2*Imax_rel*gnou_d*f_gnou_g(z3)*y_rel*c_rel^2*z4*(c_rel^2+z4^2)^(-2);

% Derivatives of concentration residuals functions
d1_k_f1_k = @(z1,z2)(1+C*(Vol*F)^(-1)*(dk_I_k1(z1) + dk_I_kr(z1) + dk_I_ks(z1,z2) + dk_I_pk(z1) + dk_I_t0(z1))*dt);

d2_na_f1_k = @(z1,z2)(C*(Vol*F)^(-1)*(-2*dna_I_nak(z2) + dna_I_ks(z1,z2))*dt);

d2_na_f2_na = @(z2)(1+C*(Vol*F)^(-1)*(dna_I_na(z2) + dna_I_bna(z2) + 3*dna_I_nak(z2) + 3*dna_I_naca(z2))*dt);

d3_ca_f2_na = (C*(Vol*F)^(-1)*3*(dca_I_naca)*dt);

d2_na_f3_ca = @(z2,z3)((-1)*C*(Vol*F)^(-1)*(dna_I_naca(z2))*dt*y_ca(z3));

d3_ca_f3_ca = @(z2,z3,z4)(1 + (C*(2*Vol*F)^(-1)*(dca_I_cal(z3) + dca_I_bca(z3) + dca_I_pca(z3) - 2*dca_I_naca) - dca_I_leak + dca_I_up(z3))*dt*y_ca(z3)...
                            + (C*(2*Vol*F)^(-1)*(I_cal(z3)+I_bca(z3)+I_pca(z3)-2*I_naca(z2,z3))-I_leak(z3,z4)+I_up(z3)-I_rel(z3,z4))*dt*c_tot*(c_buf-z3)*(c_buf+z3)*(c_tot*c_buf+(c_buf+z3)^2)^(-2));

d4_srca_f3_ca = @(z3,z4)(((-1)*dsrca_I_leak + (-1)*dsrca_I_rel(z3,z4))*dt*y_ca(z3));

d3_ca_f4_srca = @(z3,z4)(Vol*Vol_sr^(-1)*(dca_I_leak - dca_I_up(z3))*dt*y_srca(z4));

d4_srca_f4_srca = @(z3,z4)(1+Vol*Vol_sr^(-1)*((dsrca_I_leak + dsrca_I_rel(z3,z4))*y_srca(z4)+(I_leak(z3,z4)-I_up(z3)+I_rel(z3,z4))*c_srtot*(c_srbuf-z4)*(c_srbuf+z4)*(c_srtot*c_srbuf+(c_srbuf+z4)^2)^(-2))*dt);

% 4 concentrations jacobian residuals matrix
JRES = @(z1,z2,z3,z4)[d1_k_f1_k(z1,z2)  d2_na_f1_k(z1,z2) 0 0;...
                        0 d2_na_f2_na(z2) d3_ca_f2_na 0;...
                        0 d2_na_f3_ca(z2,z3) d3_ca_f3_ca(z2,z3,z4) d4_srca_f3_ca(z3,z4);...
                        0 0 d3_ca_f4_srca(z3,z4) d4_srca_f4_srca(z3,z4)];

% Initial values
vect_x = zeros(4,niter);
vect_r = zeros(1,niter);

x0 = [cant_k cant_na cant_ca cant_srca];
RES0 = RES(x0(1),x0(2),x0(3),x0(4)).';
JRES0 = JRES(x0(1),x0(2),x0(3),x0(4));

for i = 1:niter
             
    x1 = x0.' - (JRES0\RES0);
    
    vect_x(:,i) = x0.';
    vect_r(i) = abs(norm((x1.'-x0))/norm(x1.'));
    iter_fi = i;
    gnou_k1 = f_gnou_k1(x0(1));
    gnou_fca = f_gnou_fca(x0(3));
    gnou_g = f_gnou_g(x0(3));
    
    if (vect_r(i) < tol)
        break
    end
    
    x0 = x1.';
    RES0 = RES(x0(1),x0(2),x0(3),x0(4)).';
    JRES0 = JRES(x0(1),x0(2),x0(3),x0(4));
  
end

% figure(3), hold on;
% title('Error Newton iter local')
% plot (log(1:niter),log(vect_r), '*-')