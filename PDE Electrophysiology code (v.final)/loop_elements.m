function [R_V,d_V_R_V,cnou_k,cnou_na,cnou_ca,cnou_srca,...
        gnou_m,gnou_h,gnou_j,gnou_xr1,gnou_xr2,gnou_xs,...
        gnou_r,gnou_s,gnou_d,gnou_f,gnou_fca,gnou_g,gnou_k1] = loop_elements(X,R,elem,V_ant,V_nou,dt,dtt,nng,ne,nne,nnStep,I_stim,...
                                                                cant_k,cant_na,cant_ca,cant_srca,gant_m,gant_h,gant_j,gant_xr1,gant_xr2,gant_xs,...
                                                                gant_r,gant_s,gant_d,gant_f,gant_fca,gant_g)
                            

wgp     = elem.GaussWeights; 
N       = elem.N; 
Nxi     = elem.Nxi; 
Neta    = elem.Neta; 
Nzeta   = elem.Nzeta; 

R_V       = zeros(nng,1);               
d_V_R_V   = zeros(nng);

gnou_m = zeros(ne,nne);
gnou_h = zeros(ne,nne);
gnou_j = zeros(ne,nne);
gnou_d  = zeros(ne,nne);
gnou_f = zeros(ne,nne);
gnou_fca = zeros(ne,nne);
gnou_r = zeros(ne,nne);
gnou_s = zeros(ne,nne);
gnou_xs = zeros(ne,nne);
gnou_xr1 = zeros(ne,nne);
gnou_xr2 = zeros(ne,nne);
gnou_k1 = zeros(ne,nne);
gnou_g = zeros(ne,nne);

cnou_k = zeros(ne,nne);
cnou_na = zeros(ne,nne);
cnou_ca = zeros(ne,nne);
cnou_srca = zeros(ne,nne);

%Loop on elements
for ie = 1:ne
    Te = R(ie,:);
    Xe = X(Te,:);
    
    gant_m_ie = gant_m(ie,:);
    gant_h_ie = gant_h(ie,:);
    gant_j_ie = gant_j(ie,:);
    gant_d_ie = gant_d(ie,:);
    gant_f_ie = gant_f(ie,:);
    gant_fca_ie = gant_fca(ie,:);
    gant_r_ie = gant_r(ie,:);
    gant_s_ie = gant_s(ie,:);
    gant_xs_ie = gant_xs(ie,:);
    gant_xr1_ie = gant_xr1(ie,:);
    gant_xr2_ie = gant_xr2(ie,:);
    gant_g_ie = gant_g(ie,:);
    
    cant_k_ie = cant_k(ie,:);
    cant_na_ie = cant_na(ie,:);
    cant_ca_ie = cant_ca(ie,:);
    cant_srca_ie = cant_srca(ie,:);

    R_V_ie = zeros(nne,1);
    d_V_R_V_ie = zeros(nne);
    
    % Loop on element nodes
    for ine = 1:nne
        N_ine    = N(ine,:);
        Nxi_ine   = Nxi(ine,:);
        Neta_ine   = Neta(ine,:);
        Nzeta_ine   = Nzeta(ine,:);       

        Jacob = [Nxi_ine*(Xe(:,1))	Nxi_ine*(Xe(:,2))    Nxi_ine*(Xe(:,3))
                 Neta_ine*(Xe(:,1))	Neta_ine*(Xe(:,2))   Neta_ine*(Xe(:,3))
                 Nzeta_ine*(Xe(:,1))	Nzeta_ine*(Xe(:,2))  Nzeta_ine*(Xe(:,3))];
            
        dvolu = wgp(ine)*det(Jacob);
        NX_ine = Jacob\[Nxi_ine;Neta_ine; Nzeta_ine];
        
        D = 0.5;
        D_ine = D*eye(3);
       
        V_ine_ant = N_ine*V_nou(Te);
                
        gnou_m_ine = gateI_m(V_ant(R(ie,ine)),gant_m_ie(ine),dt);
        gnou_h_ine = gateI_h(V_ant(R(ie,ine)),gant_h_ie(ine),dt);
        gnou_j_ine = gateI_j(V_ant(R(ie,ine)),gant_j_ie(ine),dt);
        gnou_xr1_ine = gateI_xr1(V_ant(R(ie,ine)),gant_xr1_ie(ine),dt);
        gnou_xr2_ine = gateI_xr2(V_ant(R(ie,ine)),gant_xr2_ie(ine),dt);
        gnou_xs_ine = gateI_xs(V_ant(R(ie,ine)),gant_xs_ie(ine),dt);
        gnou_r_ine = gateI_r(V_ant(R(ie,ine)),gant_r_ie(ine),dt);
        gnou_s_ine = gateI_s(V_ant(R(ie,ine)),gant_s_ie(ine),dt);
        gnou_d_ine = gateI_d(V_ant(R(ie,ine)),gant_d_ie(ine),dt);
        gnou_f_ine = gateI_f(V_ant(R(ie,ine)),gant_f_ie(ine),dt);


        % Iterations for c_ion, gateII and I_ion variables
        [vect_x,vect_r,iter_fi,gnou_k1_ine,gnou_fca_ine,gnou_g_ine] = newton_iter_local(nnStep,10^(-9),cant_k_ie(ine),cant_na_ie(ine),cant_ca_ie(ine),cant_srca_ie(ine),...
                gnou_m_ine,gnou_h_ine,gnou_j_ine,gnou_xr1_ine,gnou_xr2_ine,gnou_xs_ine,gnou_r_ine,gnou_s_ine,gnou_d_ine,gnou_f_ine,gant_fca_ie(ine),gant_g_ie(ine),V_ant(R(ie,ine)),dtt);

        cnou_k_ine = vect_x(1,iter_fi);
        cnou_na_ine = vect_x(2,iter_fi);
        cnou_ca_ine = vect_x(3,iter_fi);
        cnou_srca_ine = vect_x(4,iter_fi);
        
        
        % Update new (nou) variables
        cnou_k(ie,ine) = cnou_k_ine;
        cnou_na(ie,ine) = cnou_na_ine;
        cnou_ca(ie,ine) = cnou_ca_ine;
        cnou_srca(ie,ine) = cnou_srca_ine;
        
        gnou_m(ie,ine) = gnou_m_ine;
        gnou_h(ie,ine) = gnou_h_ine;
        gnou_j(ie,ine) = gnou_j_ine;
        gnou_xr1(ie,ine) = gnou_xr1_ine;
        gnou_xr2(ie,ine) = gnou_xr2_ine;
        gnou_xs(ie,ine) = gnou_xs_ine;
        gnou_r(ie,ine) = gnou_r_ine;
        gnou_s(ie,ine) = gnou_s_ine;
        gnou_d(ie,ine) = gnou_d_ine;
        gnou_f(ie,ine) = gnou_f_ine;
        gnou_fca(ie,ine) = gnou_fca_ine;
        gnou_g(ie,ine) = gnou_g_ine;
        gnou_k1(ie,ine) = gnou_k1_ine;
  
        
        % Calculate f and df for the node ine of the elemet ie
        I_stim_ine = N_ine*I_stim(Te);
        [fine,dfine] = integration_point_f_df(V_ine_ant,gnou_m_ine,gnou_h_ine,gnou_j_ine,gnou_xr1_ine,gnou_xr2_ine,gnou_xs_ine,...
                    gnou_r_ine,gnou_s_ine,gnou_d_ine,gnou_f_ine,gnou_k1_ine,gnou_fca_ine,...
                    cnou_na_ine,cnou_k_ine,cnou_ca_ine,I_stim_ine);
        
   
        % Contribution node ine, elemet ie
        R_V_ie = R_V_ie + N_ine'*(N_ine*(V_nou(Te)-V_ant(Te))/dt)*dvolu...    
                      + NX_ine'*D_ine*NX_ine*V_nou(Te)*dvolu ... 
                      - N_ine'*(fine)*dvolu;
                  
        d_V_R_V_ie = d_V_R_V_ie + (N_ine'*N_ine)/dt*dvolu...
                        + (NX_ine'*D_ine*NX_ine)*dvolu...
                        - (N_ine'*N_ine)*dfine*dvolu;
                    
    end

    R_V(Te) = R_V(Te) + R_V_ie;
    d_V_R_V(Te,Te) = d_V_R_V(Te,Te) + d_V_R_V_ie;
    
end