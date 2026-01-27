% SDE: (x*(1-x) -0.5*x)dt + x\circ dW (Stratonovich)
%      (x*(1-x))dt + x*dW (Ito)

clear all
close all

warning('off','MATLAB:odearguments:RelTolIncrease');

dt = 1; % one timestep
T  = dt;

nRealizations   = 20000;
nRealizationsToPlot = 10;
plotModulus = nRealizations/nRealizationsToPlot;
nBridgeTermsMax = 25;

y0 = 0.5; % initial value

% Euler-Maruyama as 'reference'
nStepsEM = 100000;%1000000;
dtEM = dt/nStepsEM;
solEMAtT = zeros(nRealizations,1);
solEM = zeros(nStepsEM+1,1);
tEM = 0:dtEM:nStepsEM*dtEM;

% 'exact' solution by numerical integration (see 5.8 in paper draft)
nStepsEx = 100000;
dtEx = dt/nStepsEx;
tEx = 0:dtEx:nStepsEx*dtEx;
solExAtT = zeros(nRealizations,1);

for r = 1:nRealizations

    % E-M
    dW = sqrt(dtEM)*randn(1,nStepsEM);          % Brownian increments
    WEM  = [0,cumsum(dW)];                      % discretized Brownian path
%    deltaEM = WEM(end);                         % subtract linear piece to have W(T) = W(0) = 0
%    WEM_fixed = WEM - deltaEM*tEM/dt;
%    dW = [WEM_fixed(1), WEM_fixed(2:end)-WEM_fixed(1:end-1)];
    y = y0;
    solEM(1) = y;
    for j = 1:nStepsEM
        y = y + dtEM*y*(1-y) + y*dW(j);        % Ito formulation for E-M
        solEM(j+1) = y;
    end

    solEMAtT(r) = solEM(end);

    % exact
    %dWEx = sqrt(dtEx)*randn(1,nStepsEx);
    %WEx  = [0,cumsum(dWEx)];                            % discretized Brownian path
    %deltaEx = WEx(end);                        % subtract linear piece to have W(T) = W(0) = 0
    %WEx = WEx - deltaEx*tEx/dt;
    %dWEx = [WEx(1); WEx(2:end)-WEx(1:end-1)];

%    WEx = WEM_fixed;  % if number of steps is equal for exact and EM, use same random variables
    WEx = WEM;
    int_term=dtEx/2*exp(0.5)+dtEx*cumsum(exp(0.5*tEx+WEx)) - dtEx/2*exp(0.5*tEx(end)+WEx(end));
    yex = exp(0.5*tEx+WEx)./(2+int_term);
    solExAtT(r) = yex(end);

    if (mod(r,plotModulus)==0)
        figure(1);
        hold on;
        plot(tEM, solEM);
        figure(2);
        hold on;
        plot(tEx, yex);
    end
end

mean(solEMAtT)
var(solEMAtT)

mean(solExAtT)
var(solExAtT)

figure(1)
hold on
title('sample EM solutions');
xlabel('t');
ylabel('y');
figure(2)
hold on
title('sample exact solutions');
xlabel('t');
ylabel('y');


% SBB approx
xiComplete = randn(nRealizations,nBridgeTermsMax);
solAtT = zeros(nRealizations,nBridgeTermsMax);
means = zeros(nBridgeTermsMax,1);
vars  = zeros(nBridgeTermsMax,1);

options = odeset('RelTol',1e-15, 'AbsTol',1e-15);

%bridgeRange = 1:5:nBridgeTermsMax;
bridgeRange = 1:8;
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
      
      [tsol, sol] = ode45(@(t,x) SBB_ODE(t,x,dt,xi),[0 T], y0, options);
      solAtT(r,b) = sol(end);
      if (mod(r,plotModulus)==0)
        figure(3+b);
        hold on;
        plot(tsol, sol);
        title(sprintf('sample SBB ode45 solutions, b =%d', b));
        xlabel('t');
        ylabel('y');
      end
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

% figure(3+nBridgeTermsMax+1)
% plot(1:nBridgeTermsMax, means)
% figure(3+nBridgeTermsMax+2)
% plot(1:nBridgeTermsMax, vars)

figure
plot(bridgeRange, abs(means(bridgeRange)-mean(solExAtT))/abs(mean(solExAtT)))
xlabel('bridge terms');
ylabel('rel error to exact sol mean');

figure
plot(bridgeRange, abs(vars(bridgeRange)-var(solExAtT))/abs(var(solExAtT)))
xlabel('bridge terms');
ylabel('rel error to exact sol variance');


% the ODE
function  [dxdt] = SBB_ODE(t, x, dt, xi)
   sumBm =0;
   for k=2:length(xi)
      sumBm = sumBm + cos(k*pi*t/dt)*xi(k);
   end
   dbm = xi(1)/sqrt(dt) + sqrt(2/dt) .*sumBm ;
   dxdt = x*(0.5-x) + x*dbm; % Stratonovich formulation, adds -0.5x to Ito drift
end