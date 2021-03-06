clear all; close all; clc;
tic

constants

%% ---> Spatial discretization
elem = reference_element();
load('mesh9.mat')
nng = length(X(:,1)); %---> Number of global nodes (nng)
[ne, nne] = size(R); %---> Number of elements, Number of nodes per element(nne=8)(integration points)

%% ---> Time discretization
n_periodes = 1;
Time = 300;
dt = 0.5;
dtt = 0.5;


tEnd    = Time*n_periodes;     
nStep   = tEnd/dt;      %---> Number of global time steps
nnStep  = tEnd/dtt;     %---> Number of intern time steps
tt=0:dt:nStep*dt;

%% ---> Initial conditions
V0 = zeros(nng,1);
V0(:,1)= -68;
c_na = 11.6; c_k = 138.3; c_ca = 0.08*10^(-3); c_srca = 0.56;
g_m = 0; g_h = 0.75; g_j = 0.75; g_d = 0; g_f = 1; g_fca = 1; g_r = 0; g_s = 1; g_xs = 0; g_xr1 = 0; g_xr2 = 0; g_k1 = 0.5; g_g = 1;

%% ---> Iterative process
% Potential vector
V = zeros(nng,nStep); V(:,1)= V0; V_ant = V(:,1);

% Gates vectors
gant_m = g_m.*ones(ne,nne);
gant_h = g_h.*ones(ne,nne);
gant_j = g_j.*ones(ne,nne);
gant_d  = g_d.*ones(ne,nne);
gant_f = g_f.*ones(ne,nne);
gant_fca = g_fca.*ones(ne,nne);
gant_r = g_r.*ones(ne,nne);
gant_s = g_s.*ones(ne,nne);
gant_xs = g_xs.*ones(ne,nne);
gant_xr1 = g_xr1.*ones(ne,nne);
gant_xr2 = g_xr2.*ones(ne,nne);
gant_g = g_g.*ones(ne,nne);
gant_k1 = g_k1.*ones(ne,nne);

% Concentration vectors
cant_k = c_k.*ones(ne,nne);
cant_na = c_na.*ones(ne,nne);
cant_ca = c_ca.*ones(ne,nne);
cant_srca = c_srca.*ones(ne,nne);

% Save results
t = zeros (1 ,nStep/5 + 1);
V_vtk = zeros(nng, nStep/5 + 1);
cont = 2;

% time newton iterations 
for n = 2:nStep + 1
    n
    I_stim = stim(X,n,nng,dt);
    
    % loop elements
    [V_nou,cnou_k,cnou_na,cnou_ca,cnou_srca,...
    gnou_m,gnou_h,gnou_j,gnou_xr1,gnou_xr2,gnou_xs,...
    gnou_r,gnou_s,gnou_d,gnou_f,gnou_fca,gnou_g,gnou_k1] = newton_iter_time(X,R,elem,V_ant,dt,dtt,nng,ne,nne,nnStep,I_stim,...
                                                                    cant_k,cant_na,cant_ca,cant_srca,...
                                                                    gant_m,gant_h,gant_j,gant_xr1,gant_xr2,gant_xs,...
                                                                    gant_r,gant_s,gant_d,gant_f,gant_fca,gant_g);              
    V(:,n) = V_nou;
    
    % Save gates and currents to plot
    g_m_t(:,n) = gnou_m(1);
    g_h_t(:,n) = gnou_h(1);
    g_j_t(:,n) = gnou_j(1);
    g_xr1_t(:,n) = gnou_xr1(1);
    g_xr2_t(:,n) = gnou_xr2(1);
    g_xs_t(:,n) = gnou_xs(1);
    g_r_t(:,n) = gnou_r(1);
    g_s_t(:,n) = gnou_s(1);
    g_d_t(:,n) = gnou_d(1);
    g_f_t(:,n) = gnou_f(1);
    g_fca_t(:,n) = gnou_fca(1);
    g_g_t(:,n) = gnou_g(1);
    g_k1_t(:,n) = gnou_k1(1);
    g_m_t30(:,n) = gnou_m(30);
    g_h_t30(:,n) = gnou_h(30);
    g_j_t30(:,n) = gnou_j(30);
    g_xr1_t30(:,n) = gnou_xr1(30);
    g_xr2_t30(:,n) = gnou_xr2(30);
    g_xs_t30(:,n) = gnou_xs(30);
    g_r_t30(:,n) = gnou_r(30);
    g_s_t30(:,n) = gnou_s(30);
    g_d_t30(:,n) = gnou_d(30);
    g_f_t30(:,n) = gnou_f(30);
    g_fca_t30(:,n) = gnou_fca(30);
    g_g_t30(:,n) = gnou_g(30);
    g_k1_t30(:,n) = gnou_k1(30);
    
    phi_k = @(z1) phi_ion(RC,T,F,z_k,c_k0,z1);
    phi_na = @(z2) phi_ion(RC,T,F,z_na,c_na0,z2);
    phi_ca = @(z3) phi_ion(RC,T,F,z_ca,c_ca0,z3);
    phi_ks = @(z1,z2) RC*T/F*log((c_k0+p_kna*c_na0)*(z1+p_kna*z2)^(-1));
    I_na_t(n) = Cmax_na.*gnou_m(1)^3*gnou_h(1)*gnou_j(1)*(V_ant(1)-phi_na(cnou_na(1)));
    I_na_t30(n) = Cmax_na.*gnou_m(30)^3*gnou_h(30)*gnou_j(30)*(V_ant(30)-phi_na(cnou_na(30)));
    I_bna_t(n) = Cmax_bna*(V_ant(1)-phi_na(cnou_na(1)));
    I_bna_t30(n) = Cmax_bna*(V_ant(30)-phi_na(cnou_na(30)));
    I_nak_t(n) = Imax_nak*(c_k0*cnou_na(1))*((cnou_na(1)+c_nak)*(c_k0+c_kna)*(1+0.1245*exp(-0.1*V_ant(1)*F/(RC*T))+0.0353*exp(-V_ant(1)*F/(RC*T))))^(-1);
    I_nak_t30(n) = Imax_nak*(c_k0*cnou_na(30))*((cnou_na(30)+c_nak)*(c_k0+c_kna)*(1+0.1245*exp(-0.1*V_ant(30)*F/(RC*T))+0.0353*exp(-V_ant(30)*F/(RC*T))))^(-1);
    I_naca_t(n) = Imax_naca*(exp(y*V_ant(1)*F/(RC*T))*cnou_na(1)^3*c_ca0-exp((y-1)*V_ant(1)*F/(RC*T))*c_na0^3*cnou_ca(1)*y_naca)...
                *((c_naca^3+c_na0^3)*(c_cana+c_ca0)*(1+k_naca*exp((y-1)*V_ant(1)*F/(RC*T))))^(-1);
    I_naca_t30(n) = Imax_naca*(exp(y*V_ant(30)*F/(RC*T))*cnou_na(30)^3*c_ca0-exp((y-1)*V_ant(30)*F/(RC*T))*c_na0^3*cnou_ca(30)*y_naca)...
                *((c_naca^3+c_na0^3)*(c_cana+c_ca0)*(1+k_naca*exp((y-1)*V_ant(30)*F/(RC*T))))^(-1);
    I_k1_t(n) = Cmax_k1*gnou_k1(1)*(c_k0/5.4)^(1/2)*(V_ant(1)-phi_k(cnou_k(1)));
    I_k1_t30(n) = Cmax_k1*gnou_k1(30)*(c_k0/5.4)^(1/2)*(V_ant(30)-phi_k(cnou_k(30)));
    I_kr_t(n) = Cmax_kr*gnou_xr1(1)*gnou_xr2(1)*(c_k0/5.4)^(1/2)*(V_ant(1)-phi_k(cnou_k(1)));
    I_kr_t30(n) = Cmax_kr*gnou_xr1(30)*gnou_xr2(30)*(c_k0/5.4)^(1/2)*(V_ant(30)-phi_k(cnou_k(30)));
    I_ks_t(n) = Cmax_ksepi*gnou_xs(1)^2*(V_ant(1)-phi_ks(cnou_k(1),cnou_na(1)));
    I_ks_t30(n) = Cmax_ksepi*gnou_xs(30)^2*(V_ant(30)-phi_ks(cnou_k(30),cnou_na(30)));
    I_pk_t(n) =  Cmax_pk*(1+exp((25-V_ant(1))/5.98))^(-1)*(V_ant(1)-phi_k(cnou_k(1)));
    I_pk_t30(n) =  Cmax_pk*(1+exp((25-V_ant(30))/5.98))^(-1)*(V_ant(30)-phi_k(cnou_k(30)));
    I_t0_t(n) = Cmax_t0epi*gnou_r(1)*gnou_s(1)*(V_ant(1)-phi_k(cnou_k(1)));
    I_t0_t30(n) = Cmax_t0epi*gnou_r(30)*gnou_s(30)*(V_ant(30)-phi_k(cnou_k(30)));
    I_cal_t(n) = Cmax_cal*gnou_d(1)*gnou_f(1)*gnou_fca(1)*4*F^2*V_ant(1)*(RC*T)^(-1)*(cnou_ca(1)*exp(2*V_ant(1)*F*(RC*T)^(-1))-0.341*c_ca0)*(exp(2*V_ant(1)*F*(RC*T)^(-1))-1)^(-1);
    I_cal_t30(n) = Cmax_cal*gnou_d(30)*gnou_f(30)*gnou_fca(30)*4*F^2*V_ant(30)*(RC*T)^(-1)*(cnou_ca(30)*exp(2*V_ant(30)*F*(RC*T)^(-1))-0.341*c_ca0)*(exp(2*V_ant(30)*F*(RC*T)^(-1))-1)^(-1);
    I_bca_t(n) = Cmax_bca*(V_ant(1)-phi_ca(cnou_ca(1)));
    I_bca_t30(n) = Cmax_bca*(V_ant(30)-phi_ca(cnou_ca(30)));
    I_pca_t(n) = Cmax_pca*cnou_ca(1)*(c_pca+cnou_ca(1))^(-1);
    I_pca_t30(n) = Cmax_pca*cnou_ca(30)*(c_pca+cnou_ca(30))^(-1);
    I_leak_t(n) = Imax_leak*(cnou_srca(1)-cnou_ca(1));
    I_leak_t30(n) = Imax_leak*(cnou_srca(30)-cnou_ca(30));
    I_up_t(n) = Imax_up*(1+c_up^2/cnou_ca(1)^2)^(-1);
    I_up_t30(n) = Imax_up*(1+c_up^2/cnou_ca(30)^2)^(-1);
    I_rel_t(n) = Imax_rel*gnou_d(1)*gnou_g(1)*(1+y_rel*cnou_srca(1)^2*(c_rel^2+cnou_srca(1)^2)^(-1));
    I_rel_t30(n) = Imax_rel*gnou_d(30)*gnou_g(30)*(1+y_rel*cnou_srca(30)^2*(c_rel^2+cnou_srca(30)^2)^(-1));
    
    if(mod(n-1,5) < 1e-3) 
        V_vtk(:,cont) = V_nou;
        t(cont) = n*dt;
        cont = cont+1;
    end
            
  
    % Update old (ant) variables
    V_ant = V_nou;
    cant_k = cnou_k;
    cant_na = cnou_na;
    cant_ca = cnou_ca;
    cant_srca = cnou_srca;
    gant_m = gnou_m;
    gant_h = gnou_h;
    gant_j = gnou_j;
    gant_xr1 = gnou_xr1;
    gant_xr2 = gnou_xr2;
    gant_xs = gnou_xs;
    gant_r = gnou_r;
    gant_s = gnou_s;
    gant_d = gnou_d;
    gant_f = gnou_f;
    gant_fca = gnou_fca;
    gant_g = gnou_g;
    
    
end

save('pde_electrophysiology.mat', 'V_vtk', 't')
vtk_results


figure(2)
plot(tt,V(:,:))

%% ---> Plot: Gates 1
figure(3)

plot(subplot(4,4,1),tt,g_m_t(:,:))
title('g m')
xlabel('ms') 
plot(subplot(4,4,2),tt,g_h_t(:,:))
title('g h')
xlabel('ms') 
plot(subplot(4,4,3),tt,g_j_t(:,:))
title('g j')
xlabel('ms') 
plot(subplot(4,4,4),tt,g_d_t(:,:))
title('g d')
xlabel('ms')
plot(subplot(4,4,5),tt,g_f_t(:,:))
title('g f')
xlabel('ms') 
plot(subplot(4,4,6),tt,g_fca_t(:,:))
title('g fca')
xlabel('ms') 
plot(subplot(4,4,7),tt,g_r_t(:,:))
title('g r')
xlabel('ms') 
plot(subplot(4,4,8),tt,g_s_t(:,:))
title('g s')
xlabel('ms') 
plot(subplot(4,4,9),tt,g_xs_t(:,:))
title('g xs')
xlabel('ms') 
plot(subplot(4,4,10),tt,g_xr1_t(:,:))
title('g xr1')
xlabel('ms') 
plot(subplot(4,4,11),tt,g_xr2_t(:,:))
title('g xr2')
xlabel('ms') 
plot(subplot(4,4,12),tt,g_k1_t(:,:))
title('g k1')
xlabel('ms') 
plot(subplot(4,4,13),tt,g_g_t(:,:))
title('g g')
xlabel('ms') 

%% ---> Plot: Gates 30
figure(4)

plot(subplot(4,4,1),tt,g_m_t30(:,:))
title('g m')
xlabel('ms') 
plot(subplot(4,4,2),tt,g_h_t30(:,:))
title('g h')
xlabel('ms') 
plot(subplot(4,4,3),tt,g_j_t30(:,:))
title('g j')
xlabel('ms') 
plot(subplot(4,4,4),tt,g_d_t30(:,:))
title('g d')
xlabel('ms')
plot(subplot(4,4,5),tt,g_f_t30(:,:))
title('g f')
xlabel('ms') 
plot(subplot(4,4,6),tt,g_fca_t30(:,:))
title('g fca')
xlabel('ms') 
plot(subplot(4,4,7),tt,g_r_t30(:,:))
title('g r')
xlabel('ms') 
plot(subplot(4,4,8),tt,g_s_t30(:,:))
title('g s')
xlabel('ms') 
plot(subplot(4,4,9),tt,g_xs_t30(:,:))
title('g xs')
xlabel('ms') 
plot(subplot(4,4,10),tt,g_xr1_t30(:,:))
title('g xr1')
xlabel('ms') 
plot(subplot(4,4,11),tt,g_xr2_t30(:,:))
title('g xr2')
xlabel('ms') 
plot(subplot(4,4,12),tt,g_k1_t30(:,:))
title('g k1')
xlabel('ms') 
plot(subplot(4,4,13),tt,g_g_t30(:,:))
title('g g')
xlabel('ms')

%% ---> Plot: Currents 1
figure(5)

plot(subplot(4,4,1),tt,I_na_t(:))
title('I na')
xlabel('ms') 
ylabel('pA/pF') 

plot(subplot(4,4,2),tt,I_bna_t(:))
title('I bna')
xlabel('ms') 
ylabel('pA/pF')

plot(subplot(4,4,3),tt,I_nak_t(:))
title('I nak')
xlabel('ms') 
ylabel('pA/pF')

plot(subplot(4,4,4),tt,I_naca_t(:))
title('I naca')
xlabel('ms') 
ylabel('pA/pF')

plot(subplot(4,4,5),tt,I_k1_t(:))
title('I k1')
xlabel('ms') 
ylabel('pA/pF')

plot(subplot(4,4,6),tt,I_kr_t(:))
title('I kr')
xlabel('ms') 
ylabel('pA/pF')

plot(subplot(4,4,7),tt,I_ks_t(:))
title('I ks')
xlabel('ms') 
ylabel('pA/pF')

plot(subplot(4,4,8),tt,I_pk_t(:))
title('I pk')
xlabel('ms') 
ylabel('pA/pF')


plot(subplot(4,4,9),tt,I_t0_t(:))
title('I to')
xlabel('ms') 
ylabel('pA/pF')

plot(subplot(4,4,10),tt,I_cal_t(:))
title('I cal')
xlabel('ms') 
ylabel('pA/pF')

plot(subplot(4,4,11),tt,I_bca_t(:))
title('I bca')
xlabel('ms') 
ylabel('pA/pF')

plot(subplot(4,4,12),tt,I_pca_t(:))
title('I pca')
xlabel('ms') 
ylabel('pA/pF')

plot(subplot(4,4,13),tt,I_leak_t(:))
title('I leak')
xlabel('ms') 
ylabel('mM/ms')

plot(subplot(4,4,14),tt,I_up_t(:))
title('I up')
xlabel('ms') 
ylabel('mM/ms')

plot(subplot(4,4,15),tt,I_rel_t(:))
title('I rel')
xlabel('ms') 
ylabel('mM/ms')


%% ---> Plot: Currents 30
figure(6)

plot(subplot(4,4,1),tt,I_na_t30(:))
title('I na')
xlabel('ms') 
ylabel('pA/pF') 

plot(subplot(4,4,2),tt,I_bna_t30(:))
title('I bna')
xlabel('ms') 
ylabel('pA/pF')

plot(subplot(4,4,3),tt,I_nak_t30(:))
title('I nak')
xlabel('ms') 
ylabel('pA/pF')

plot(subplot(4,4,4),tt,I_naca_t30(:))
title('I naca')
xlabel('ms') 
ylabel('pA/pF')

plot(subplot(4,4,5),tt,I_k1_t30(:))
title('I k1')
xlabel('ms') 
ylabel('pA/pF')

plot(subplot(4,4,6),tt,I_kr_t30(:))
title('I kr')
xlabel('ms') 
ylabel('pA/pF')

plot(subplot(4,4,7),tt,I_ks_t30(:))
title('I ks')
xlabel('ms') 
ylabel('pA/pF')

plot(subplot(4,4,8),tt,I_pk_t30(:))
title('I pk')
xlabel('ms') 
ylabel('pA/pF')


plot(subplot(4,4,9),tt,I_t0_t30(:))
title('I to')
xlabel('ms') 
ylabel('pA/pF')

plot(subplot(4,4,10),tt,I_cal_t30(:))
title('I cal')
xlabel('ms') 
ylabel('pA/pF')

plot(subplot(4,4,11),tt,I_bca_t30(:))
title('I bca')
xlabel('ms') 
ylabel('pA/pF')

plot(subplot(4,4,12),tt,I_pca_t30(:))
title('I pca')
xlabel('ms') 
ylabel('pA/pF')

plot(subplot(4,4,13),tt,I_leak_t30(:))
title('I leak')
xlabel('ms') 
ylabel('mM/ms')

plot(subplot(4,4,14),tt,I_up_t30(:))
title('I up')
xlabel('ms') 
ylabel('mM/ms')

plot(subplot(4,4,15),tt,I_rel_t30(:))
title('I rel')
xlabel('ms') 
ylabel('mM/ms')


