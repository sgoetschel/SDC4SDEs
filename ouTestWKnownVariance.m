% SDE: kappa(mu-x)dt + sigma \circ dW (Stratonovich)
%      kappa(mu-x)dt + sigma dW (Ito)

clear all
close all

seed = 2348;
rng(seed)

warning('off','MATLAB:odearguments:RelTolIncrease');

dt = 2; % one timestep
T  = dt;

nRealizations   = 10000;
nRealizationsToPlot = 10;
plotModulus = nRealizations/nRealizationsToPlot;
nBridgeTermsMax = 2;

y0 = 0.5; % initial value
kappa = 0.5;
mu=0;
sigma=1;

exactMean = mu-(mu-y0)*exp(-kappa*T)
exactVar  = sigma^2/(2*kappa)*(1-exp(-2*kappa*T))

% SBB approx
xiComplete = randn(nRealizations,nBridgeTermsMax);
solAtT = zeros(nRealizations,nBridgeTermsMax);
means = zeros(nBridgeTermsMax,1);
vars  = zeros(nBridgeTermsMax,1);

options = odeset('RelTol',1e-15, 'AbsTol',1e-15);

%bridgeRange = 1:5:nBridgeTermsMax;
bridgeRange = 1:2;
for b = bridgeRange
   for r = 1:nRealizations
      xi = xiComplete(r,1:b);
      % how to fix the endpoint of the Wiener process here? 
      % just set eta_0 = xi(1) to zero.
%      xi(1) = 0;
      
%       time = 0:0.000001:T;
%       sumBm =0;
%       for k=1:length(xi)-1
%          sumBm = sumBm + sin(k*pi*time/dt)*xi(k+1)/k;
%       end
%       bm = time*xi(1)/sqrt(dt) + sqrt(2*dt)/pi .*sumBm;
% 
%       figure(30)
%       hold on
%       plot(time, bm)
      
      [tsol, sol] = ode45(@(t,x) SBB_ODE(t,x,dt,xi,kappa,mu,sigma),[0 T], y0, options);
      solAtT(r,b) = sol(end);
%       if (mod(r,plotModulus)==0)
%         figure(3+b);
%         hold on;
%         plot(tsol, sol);
%         title(sprintf('sample SBB ode45 solutions, b =%d', b));
%         xlabel('t');
%         ylabel('y');
%       end
   end
   mb = mean(solAtT(:,b))
   vb = var(solAtT(:,b))
   means(b) = mean(solAtT(:,b));
   vars(b) = var(solAtT(:,b));
end

% what does the variance of the solution mean anyway? 

% observation: starting from 5 bridge terms, ode45 gives the same solution
% regardless whether more bridge terms are used? WHY? temporal resolution
% not fine enough to resolve oszillations anymore; with tolerance 1e-15 in
% ode solver this gives slight differences. The error goes down in mean and
% variance, but only a little.

figure
plot(bridgeRange, abs(means(bridgeRange)-exactMean),'o')
xlabel('bridge terms');
ylabel('rel error to exact sol mean');

figure
plot(bridgeRange, abs(vars(bridgeRange)-exactVar),'x')
xlabel('bridge terms');
ylabel('rel error to exact sol variance');


% the ODE
function  [dxdt] = SBB_ODE(t, x, dt, xi, kappa, mu, sigma)
   sumBm =0;
   for k=2:length(xi)
      sumBm = sumBm + cos(k*pi*t/dt)*xi(k);
   end
   dbm = xi(1)/sqrt(dt) + sqrt(2/dt) .*sumBm ;
   dxdt = kappa*(mu-x) + sigma*dbm; 
end