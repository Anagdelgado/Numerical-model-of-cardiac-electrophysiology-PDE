function [I_stim] = stim(X,n,nng,dt)
I_stim = zeros(nng,1);
for i = 1:nng
    if (abs(X(i,1)) < 1e-5) && n*dt < 5
        I_stim(i) = 20;
    else
        I_stim(i)=0;
    end
end