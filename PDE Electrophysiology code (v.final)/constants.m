% Initial concentrations
c_na0 = 140; c_k0 = 5.4; c_ca0 = 2;

%Elementary charge per ion
z_na = 1; z_k = 1; z_ca = 2;

% Maximum currents (in pA/pF=volt/s, o mV/ms)
Imax_naca = 1000.; Imax_nak = 1.362; 

%(in mM/s)
Imax_up = 0.425/1000; Imax_rel = 8.232/1000; Imax_leak = 0.08/1000.; 

% Maximum conductances (in nS/pF)
Cmax_na = 14.838; Cmax_bna = 0.00029; Cmax_bca = 0.000592; Cmax_pca = 0.825;
Cmax_k1 = 5.405; Cmax_kr = 0.096; %according to Tusscher et al
Cmax_ksepi = 0.245; Cmax_ksendo = 0.245; Cmax_ksm = 0.062;
Cmax_pk = 0.0146;  Cmax_t0epi = 0.294;  Cmax_t0endo = 0.073; Cmax_t0m = 0.294;

%(in mm3/microF s)
Cmax_cal = 0.175;

% Half saturation constants
c_cana = 1.38; c_naca = 87.5; c_kna = 1; c_nak = 40;
c_pca = 0.0005; c_up = 0.00025; c_rel = 0.25; c_buf = 0.001; c_srbuf = 0.3;

% Other parameters
k_naca = 0.1; p_kna = 0.03; c_tot = 0.15; c_srtot = 10;
y_naca = 2.5; y = 0.35; y_rel = 2;

% Generic constants
RC = 8.3143; F = 96.4867; T = 310;
C = 185; Vol = 16404; Vol_sr = 1094;  


