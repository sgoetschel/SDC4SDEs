% SDE: kappa(mu-x)dt + sigma \circ dW (Stratonovich)
%      kappa(mu-x)dt + sigma dW (Ito)

clear all
close all

format long

seed = 2348;
rng(seed)

T = 1/8; 

dts = [2; 1; 1/2; 1/4; 1/8; ];% 1/16; 1/32];% 1/256];

NNfinest = T/dts(end) % 1; %one timestep
dtFinest  = T/NNfinest; % one timestep

nDifferentSteps = length(dts);

nBatches = 1;
nRealizations   = 100000;
nRealizationsToPlot = 10;
plotModulus = nRealizations/nRealizationsToPlot;
%nBridgeTermsMax = 25;

order = 1; % number of equations/unknowns

y0 = 0.5; % initial value
kappa = 0.5;
mu=0;
sigma=1;

exactMean = mu-(mu-y0)*exp(-kappa*T)
exactVar  = sigma^2/(2*kappa)*(1-exp(-2*kappa*T))
% comparing to these exact quantities we additionally have a sample error.
% not just the time/SBB approximation errors
% use E-M with massively many time steps as a bechmark solution? Weak order
% 1 to approximate the mean

% Euler-Maruyama as 'reference'
nStepsEM = 128000;
dtEM = T/nStepsEM;
solEMAtT = zeros(nBatches,nRealizations);

meansEM = zeros(nDifferentSteps, nBatches);
varsEM = zeros(nDifferentSteps, nBatches);
varEM  = zeros(nDifferentSteps, 1);
meanEM = zeros(nDifferentSteps, 1);

% SBB approx

bridgeRange = 1:2;
nBridgeTermsMax = bridgeRange(end);

nRndVar = NNfinest;
eta = cell(1,order);
xi = cell(1,order);
% random variable for linear part of bridge
for k=1:order
    eta{k} = cumsum(randn(nRealizations,nRndVar),2).*sqrt(T)./sqrt(nRndVar); 
end

%draw samples for more than one Karhunen-Loeve expansion term
for k=1:order
    xi{k} = randn(nRealizations,nBridgeTermsMax-1, NNfinest);
end
% for k=1:order
%     eta{k} = cumsum(xi{k}(:,1,:),3).*sqrt(dtFinest);%cumsum(randn(nRealizations,nRndVar),2).*sqrt(T)./sqrt(nRndVar); 
% end


%%%% Analytical solution of RODE at T%%%%
bridgeRangeAna =[1,2,5,10,100,1000,10000];
nBridgeTermsAna = 100;
solAnaAtT = zeros(nRealizations,1);
meanAna = zeros(length(bridgeRangeAna),1);
varAna = zeros(length(bridgeRangeAna),1);
t=T;
DeltaT = T;
counter = 1;
for nBridgeTermsAna = bridgeRangeAna
for l=1:nRealizations
    eta0 = randn(1,1);
    eta = randn(nBridgeTermsAna-1,1);
    a = -kappa;
    b = kappa * mu + sigma * eta0/sqrt(DeltaT);
    c = sigma * sqrt(2/DeltaT);
    sumForIC = 0;
    sumBridgeTerms = 0;
    for j=1:nBridgeTermsAna-1
        omega = j*pi/DeltaT;
        sumForIC = sumForIC + eta(j)/(a^2+omega^2);
        term = eta(j)*exp(-a*t)*(omega*sin(omega*t) - a*cos(omega*t))/(a^2+omega^2);
        sumBridgeTerms = sumBridgeTerms + term;
    end
    K = y0 + a*c*sumForIC+b/a;
    solAnaAtT(l) = exp(a*t)*(c*sumBridgeTerms-b*exp(-a*t)/a + K);
end
meanAna(counter)=mean(solAnaAtT);
varAna(counter)=var(solAnaAtT);
counter=counter+1;
end

for i=1:nDifferentSteps
    dtEM=dts(i)
    nStepsEM = T/dtEM
    for batchNo = 1:nBatches
        for l=1:nRealizations
            solEM = zeros(nStepsEM+1,1);
            %WEM = [0 eta{1}(l,:)];
            %dW = WEM(2:end)-WEM(1:end-1);
            dW = randn(1,nStepsEM).*sqrt(T)./sqrt(nStepsEM);          % Brownian increments
            WEM  = [0,cumsum(dW)];      
            y = y0;
            solEM(1) = y0;
            for j = 1:nStepsEM
                y = y + dtEM*kappa*(mu-y) + sigma*dW(j);        % Ito formulation for E-M
                solEM(j+1) = y;
            end
            solEMAtT(batchNo,l) = y;
            
%             if (mod(l,plotModulus)==0)
%                 figure(i);
%                  hold on;
%                 plot(0:dtEM:T, solEM);
%             end
        end
        meansEM(i,batchNo) = mean(solEMAtT(batchNo,:));
        varsEM(i,batchNo)  = var(solEMAtT(batchNo,:));
    end
    meanEM(i) = sum(meansEM(i,:))/nBatches;
    varEM(i)  = var(solEMAtT(:));
end

figure;
hold on;
semilogy(dts, abs(meanEM(:)-exactMean));

%varEM=var(solEMAtT)
%meanEM=mean(solEMAtT)

%xiComplete = randn(nRealizations,nBridgeTermsMax);

means = zeros(nBridgeTermsMax,nDifferentSteps);
vars  = zeros(nBridgeTermsMax,nDifferentSteps);

errMean = zeros(nBridgeTermsMax, nDifferentSteps);
errVar  = zeros(nBridgeTermsMax, nDifferentSteps);

% SDC setup
tol=1e-12;
col_points=4;
max_sweeps=100;
nodes = 'lobatto';
strInit = 'eulerSBBInit';

for i=1:nDifferentSteps
    
dt=dts(i)

solsAtT = zeros(nRealizations,nBridgeTermsMax);
for b = bridgeRange
   [solsAtT(:,b), ~] = sdeSolverSBBSDC(tol, NNfinest, dt, 0, T, col_points, max_sweeps, nRealizations, y0, nodes, sigma, ...
       @(y)drift_Stra(y,kappa,mu), {@(y)diffusion(y,sigma)}, @(y)drift_Ito(y,kappa,mu), order, eta, xi, strInit, b);
   mb = mean(solsAtT(:,b))
   vb = var(solsAtT(:,b))
   means(b,i) = mb;
   vars(b,i) = vb;
end

%errMean(:,i) = abs(means(:,i)-exactMean)/abs(exactMean);
%errVar(:,i) =  abs(vars(:,i)-exactVar)/abs(exactVar);

% errMean(:,i) = abs(means(:,i)-meanEM)/abs(meanEM);
% errVar(:,i) =  abs(vars(:,i)-varEM)/abs(varEM);

% figure
% plot(bridgeRange, errMean(:,i))
% xlabel('bridge terms');
% ylabel('rel error to exact sol mean');
% title(['dt = ', num2str(dt)]);
% 
% 
% figure
% plot(bridgeRange, errVar(:,i))
% xlabel('bridge terms');
% ylabel('rel error to exact sol variance');
% title(['dt = ', num2str(dt)]);

end

plot(dts, abs(means(1,:)-exactMean));
plot(dts, abs(means(2,:)-exactMean));
legend('EM','SBB-SDC-4 m=1', 'SBB-SDC-4 m=2')
% 
% figure
% hold on
% for b = bridgeRange
%     plot( dts, errMean(b,:), 'DisplayName', ['b = ', num2str(b)]);
% end
% xlabel('dt');
% ylabel('rel error to exact sol mean');
% legend
% 
% figure
% hold on
% for b = bridgeRange
%     plot( dts, errVar(b,:), 'DisplayName', ['b = ', num2str(b)]);
% end
% xlabel('dt');
% ylabel('rel error to exact sol variance');
% legend
% 
% % what does the variance of the solution mean anyway? 
% 
% % observation: starting from 5 bridge terms, ode45 gives the same solution
% % regardless whether more bridge terms are used? WHY? temporal resolution
% % not fine enough to resolve oszillations anymore; with tolerance 1e-15 in
% % ode solver this gives slight differences. The error goes down in mean and
% % variance, but only a little.
% 
% % figure
% % plot(bridgeRange, abs(means(bridgeRange)-exactMean)/abs(exactMean))
% % xlabel('bridge terms');
% % ylabel('rel error to exact sol mean');
% % 
% % figure
% % plot(bridgeRange, abs(vars(bridgeRange)-exactVar)/abs(exactVar))
% % xlabel('bridge terms');
% % ylabel('rel error to exact sol variance');
% 

function [r] = drift_Stra(y,kappa,mu)
    r = kappa*(mu-y);
end

function [r] = drift_Ito(y,kappa,mu)
    r = kappa*(mu-y);
end

function [r] = diffusion(y,sigma)
    r = sigma*ones(1,length(y));
end

% % the ODE
% function  [dxdt] = SBB_ODE(t, x, dt, xi, kappa, mu, sigma)
%    sumBm =0;
%    for k=2:length(xi)
%       sumBm = sumBm + cos(k*pi*t/dt)*xi(k);
%    end
%    dbm = xi(1)/sqrt(dt) + sqrt(2/dt) .*sumBm ;
%    dxdt = kappa*(mu-x) + sigma*dbm; 
% end