function element = reference_element()

Xe_ref  = [-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1];
zgp     = (1/sqrt(3))*[-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1];
wgp     = ones(1,8);
ngaus   = length(wgp);
xi      = zgp(:,1);
eta     = zgp(:,2);
zeta    = zgp(:,3);
N       =  (1/8)*[(1-xi).*(1-eta).*(1-zeta), (1+xi).*(1-eta).*(1-zeta), (1+xi).*(1+eta).*(1-zeta), (1-xi).*(1+eta).*(1-zeta), (1-xi).*(1-eta).*(1+zeta), (1+xi).*(1-eta).*(1+zeta), (1+xi).*(1+eta).*(1+zeta), (1-xi).*(1+eta).*(1+zeta)];
Nxi     =  (1/8)*[-(1-eta).*(1-zeta), (1-eta).*(1-zeta), (1+eta).*(1-zeta), -(1+eta).*(1-zeta), -(1-eta).*(1+zeta), (1-eta).*(1+zeta), (1+eta).*(1+zeta), -(1+eta).*(1+zeta)];
Neta    =  (1/8)*[-(1-xi).*(1-zeta), -(1+xi).*(1-zeta), (1+xi).*(1-zeta), (1-xi).*(1-zeta), -(1-xi).*(1+zeta), -(1+xi).*(1+zeta), (1+xi).*(1+zeta), (1-xi).*(1+zeta)];
Nzeta   =  (1/8)*[-(1-xi).*(1-eta), -(1+xi).*(1-eta), -(1+xi).*(1+eta), -(1-xi).*(1+eta), (1-xi).*(1-eta), (1+xi).*(1-eta), (1+xi).*(1+eta), (1-xi).*(1+eta)];

element.degree          = 1;
element.Xe_ref          = Xe_ref; 
element.ngaus           = ngaus; 
element.GaussPoints     = zgp;
element.GaussWeights    = wgp;
element.N               = N; 
element.Nxi             = Nxi; 
element.Neta            = Neta; 
element.Nzeta           = Nzeta;  

end