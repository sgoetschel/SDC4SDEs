function  [dxdt] = SBB_ODE_exp(s,x, eta, deltaT, xi, lambda, beta, thisTimeStep)

sumBm =0;
for k=1:length(xi)
    sumBm = sumBm + cos(k*pi*(s)/deltaT)*xi(k);%fix this for multidim sde
end
dbm = eta/sqrt(deltaT) + sqrt(2/deltaT) .*sumBm;

dxdt = (lambda-0.5*beta^2)*x +  beta*x*dbm;

end