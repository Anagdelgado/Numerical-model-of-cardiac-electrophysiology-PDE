function [V_nou,cnou_k,cnou_na,cnou_ca,cnou_srca,...
    gnou_m,gnou_h,gnou_j,gnou_xr1,gnou_xr2,gnou_xs,...
    gnou_r,gnou_s,gnou_d,gnou_f,gnou_fca,gnou_g,gnou_k1] = newton_iter_time(X,R,elem,V_ant,dt,dtt,nng,ne,nne,nnStep,I_stim,...
                                                                            cant_k,cant_na,cant_ca,cant_srca,...
                                                                            gant_m,gant_h,gant_j,gant_xr1,gant_xr2,gant_xs,...
                                                                            gant_r,gant_s,gant_d,gant_f,gant_fca,gant_g)                                                                   
tol = 1e-10;
cont = 1;
V_nou = V_ant;
error_plot_V =[];
error_plot_dV = [];
V_error = 1;

while (cont < 400) && ((V_error > tol) || (dV_error > tol))                                                         
    [R_V,d_V_R_V,cnou_k,cnou_na,cnou_ca,cnou_srca,...
    gnou_m,gnou_h,gnou_j,gnou_xr1,gnou_xr2,gnou_xs,...
    gnou_r,gnou_s,gnou_d,gnou_f,gnou_fca,gnou_g,gnou_k1] = loop_elements(X,R,elem,V_ant,V_nou,dt,dtt,nng,ne,nne,nnStep,I_stim,...
                                                                cant_k,cant_na,cant_ca,cant_srca,gant_m,gant_h,gant_j,gant_xr1,gant_xr2,gant_xs,...
                                                              gant_r,gant_s,gant_d,gant_f,gant_fca,gant_g);
    
                                                      
    dV = - d_V_R_V\R_V;
    V_nou = V_nou + dV;
    V_error = norm(dV)/norm(V_nou);
    dV_error = norm(R_V);
    error_plot_V =[error_plot_V V_error];
    error_plot_dV =[error_plot_dV dV_error];
    cont = cont +1;
    
% Confirmar que no se deben actualizar
%     cant_k = cnou_k;
%     cant_na = cnou_na;
%     cant_ca = cnou_ca;
%     cant_srca = cnou_srca;
%     gant_m = gnou_m;
%     gant_h = gnou_h;
%     gant_j = gnou_j;
%     gant_xr1 = gnou_xr1;
%     gant_xr2 = gnou_xr2;
%     gant_xs = gnou_xs;
%     gant_r = gnou_r;
%     gant_s = gnou_s;
%     gant_d = gnou_d;
%     gant_f = gnou_f;
%     gant_fca = gnou_fca;
%     gant_g = gnou_g;
end

semilogy([1:cont-1],error_plot_V,'ro-',[1:cont-1],error_plot_dV,'bo-');

% figure(1), hold on;
% title('Error Newton iter global')
% plot (log(1:cont-1),log(error_plot_V), '*-')


if (V_error > tol) 
    fprintf('Newton no convergeix \n')
    return
end

end
                                                        
