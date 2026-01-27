function  [dxdt] = SBB_ODE_TP3(s,x, eta, deltaT, xi, lambda, beta)

sumBm =0;
for k=1:length(xi)
    sumBm = sumBm + cos(k*pi*s/deltaT)*xi(k);%fix this for multidim sde
end
dbm = eta/sqrt(deltaT) + sqrt(2/deltaT) .*sumBm;

dxdt = lambda*x*(0.5-x) + beta*x*dbm;

end