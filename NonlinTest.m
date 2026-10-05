%% ================================================================
%  Convergence study for a nonlinear SDE with Wong-Zakai approximation
%
%  SDE:
%
%      dX = lambda*x.*(1-x) dt + beta x dW              (Ito)
%      dX = lambda*x.*(1-x)-0.5*beta*x dt + beta x dW   (Stratonovich)
%
%  Three levels:
%
%      X_h^K : RK4 approximation of Wong-Zakai ODE
%      X^K   : exact solution of Wong-Zakai ODE
%      X     : exact OU solution
%
%  Error decomposition:
%
%      X_h^K - X = (X_h^K-X^K) + (X^K-X)
%
%  Three plots:
%
%      1. RK4 time-step convergence
%      2. Wong-Zakai convergence
%      3. Total error
%
% ================================================================

clearvars;
close all;
clc;

%% Parameters
T      = 1;
X0     = 0.5;
lambda = 1;
beta   = 1;

nStepsRef     = 1024;                    % max number of time steps at which Brownian motion is sampled
nSteps        = [4, 8, 16, 32, 64, 128]; % time steps for convergence study
nStepsizes    = length(nSteps);
nRealizations = 10000;

tRef = linspace(0,T,nStepsRef+1);
dtRef = tRef(2)-tRef(1);

rng(12345);

%% ================================================================
% Fourier/Wong-Zakai representation
% ================================================================
%
% We use
%
%   dW^K/dt =
%       xi_0/sqrt(T)
%       + sqrt(2/T) sum_{k=1}^K xi_k cos(k*pi*t/T)
%
% with xi_k ~ N(0,1).

Kvalues = [1 2 4 8 16 32 64];
Kmax = max(Kvalues);
nK = length(Kvalues);


%% ================================================================
% Analytical solution
% ================================================================

solRef = zeros(nRealizations, nStepsRef+1);
solW_K = zeros(nRealizations, nK, nStepsRef+1);

% sample Brownian motion at finest resolution in time
eta = zeros(nRealizations, nStepsRef+1);
eta(:, 2:end) = cumsum(randn(nRealizations,nStepsRef),2).*sqrt(T)./sqrt(nStepsRef);

dW = eta(:,2:end)-eta(:,1:end-1);  % increments of Brownian motion
xiEndpoint = eta(:,end);           % end value of Brownian motion after nStepsRef time steps

% compute Wong-Zakai approximation *of the Brownian Bridge connecting 0 to xiEndpoint=W(T)*
k = 1:Kmax-1;                       % without the linear part
C = cos(pi/T * tRef(1:end-1)' * k); %cos(pi/nStepsRef * n * k);
xi = (sqrt(2/T) * C' * dW')';
% Reconstruction matrix
S = sin(pi/T * tRef' * k);
% Brownian-bridge coefficients
a = sqrt(2*T) ./ (pi*k) .* xi;

%% ================================================================
%  Evaluate Wong-Zakai approximation error for the full solution
%  No use of Brownian bridges here!
%  Reference: with sampled Brownian motion on fine time grid
%  Approximation: use KL-expansion instead of sampled Brownian motion
for i = 1:nK
    K = Kvalues(i);

    for  l=1:nRealizations
        W_K = (tRef/T).*xiEndpoint(l,:) + (S(:,1:K-1)*a(l,1:K-1)')'; 

        % plot comparing one realization of Brownian motion and its BB-Fourier
        % approximation
        % figure;
        % plot(tRef, eta(l,:), tRef, W_K) 
    
        solRef(l,:) = exactTP3(lambda, beta, 0, T, dtRef, X0, eta(l,:));
        solW_K(l,i,:) = exactTP3(lambda, beta, 0, T, dtRef, X0, W_K);
    end
end

meanRef = mean(solRef(:,end));
varRef  = var(solRef(:,end));

fprintf('\nReference statistics:\n');
fprintf('  E[X(T)]     = %.12e\n',meanRef);
fprintf('  Var[X(T)]   = %.12e\n',varRef);

for i = 1:nK
    K = Kvalues(i);
    fprintf('\n----------------------------------------\n');
    fprintf('Number of Fourier modes: %d\n',K);
    
    meanApprox(i) = mean(solW_K(:,i,end));
    varApprox(i)  = var(solW_K(:,i,end));
    fprintf('  Mean             = %.12e\n',meanApprox(i));
    fprintf('  Variance         = %.12e\n',varApprox(i));

    % strong final time errors
    err = solW_K(:,i, end) - solRef(:,end);

    errorWZ_L1(i) = mean(abs(err));
    errorWZ_L2(i) = sqrt(mean(err.^2));

    % weak final time errors
    weakMean(i) = abs(meanApprox(i)-meanRef);
    weakVar(i) = abs(varApprox(i)-varRef);

    fprintf('Strong L1 error  = %.12e\n',errorWZ_L1(i));
    fprintf('Strong L2 error  = %.12e\n',errorWZ_L2(i));

    fprintf('Weak mean error  = %.12e\n',weakMean(i));
    fprintf('Weak var error   = %.12e\n',weakVar(i));
end

% Fit slopes with respect to K.
%
% If error ~ K^(-q), then
%
%   log(error) = -q log(K) + C.

pWZ_L1 = polyfit(log(Kvalues),log(errorWZ_L1'),1);
pWZ_L2 = polyfit(log(Kvalues),log(errorWZ_L2'),1);

fprintf('\n');
fprintf('====================================================\n');
fprintf('WONG-ZAKAI CONVERGENCE\n');
fprintf('====================================================\n');
fprintf('L1 order with respect to K = %.4f\n',-pWZ_L1(1));
fprintf('L2 order with respect to K = %.4f\n',-pWZ_L2(1));

%% Plot Wong-Zakai convergence

figure;

loglog(Kvalues,errorWZ_L1,'o-','DisplayName','L^1 error');
hold on;
loglog(Kvalues,errorWZ_L2,'s-','DisplayName','L^2 error');

grid on;

xlabel('number of Fourier modes K');
ylabel('error at T');

title('Wong-Zakai convergence');

legend('Location','best');



%% ================================================================
% RK4 TIME-STEP CONVERGENCE
% ================================================================
%
% Fix K and vary h.
%
% The reference is X^K, which is known analytically.
%
% Therefore this measures ONLY:
%
%       X_h^K - X^K

Kfixed = 64;

hValues = T/(Kfixed) ./ 2.^(0:6);

nH = length(hValues);

errorTimeL1 = zeros(nH,1);
errorTimeL2 = zeros(nH,1);

idxK = find(Kvalues == Kfixed);

for ih = 1:nH

    h = hValues(ih);

    Xh = zeros(nRealizations,1);

    for r = 1:nRealizations

        Xh(r) = solve_WZ_RK4( ...
            X0,T,h,lambda,beta,...
            xiEndpoint(r), xi(r,1:Kfixed-1),Kfixed);
    end

    err_to_Ref = Xh - solRef(:,end);
    err  = Xh - solW_K(:,idxK,end);

    errorTimeL1(ih) = mean(abs(err));

    errorTimeL2(ih) = sqrt(mean(err.^2));

    errorRefTimeL1(ih) = mean(abs(err_to_Ref));
    errorRefTimeL2(ih) = sqrt(mean(err_to_Ref.^2));
end

% Estimate slopes
pTimeL1 = polyfit(log(hValues),log(errorTimeL1'),1);
pTimeL2 = polyfit(log(hValues),log(errorTimeL2'),1);

fprintf('\n');
fprintf('====================================================\n');
fprintf('RK4 TIME-STEP CONVERGENCE\n');
fprintf('====================================================\n');
fprintf('K = %d\n',Kfixed);
fprintf('L1 order = %.4f\n',pTimeL1(1));
fprintf('L2 order = %.4f\n',pTimeL2(1));

%% Plot 1: time-step convergence

figure;

loglog(hValues,errorTimeL1,'o-','DisplayName','L^1 error to WZ');
hold on;
loglog(hValues,errorTimeL2,'s-','DisplayName','L^2 error to WZ');
loglog(hValues,errorRefTimeL1,'s-','DisplayName','L^1 error to ref');
loglog(hValues,errorRefTimeL2,'s-','DisplayName','L^2 error to ref');

grid on;

xlabel('time step h');
ylabel('error at T');

title(sprintf('RK4 time-step convergence, K = %d',Kfixed));

legend('Location','best');

%% ================================================================
% 3. TOTAL ERROR
% ================================================================
%
% Now solve the Wong-Zakai ODE numerically:
%
%       X_h^K
%
% and compare directly against the exact OU solution:
%
%       X.
%
% Thus:
%
%       total error = X_h^K - X.
%
% We choose h proportional to 1/K and sufficiently small such
% that the time discretization error is small compared with the
% Wong-Zakai error.

C = 80;

errorTotal_L2 = zeros(nK,1);
errorTime_L2 = zeros(nK,1);

for i = 1:nK

    K = Kvalues(i);

    % Resolve the highest Fourier frequency increasingly well
    h = T/(C*K);

    Xh = zeros(nRealizations,1);

    for r = 1:nRealizations

        Xh(r) = solve_WZ_RK4( ...
            X0,T,h,lambda,beta,...
            xiEndpoint(r), xi(r,1:K-1),K);

    end

    % Time discretization component
    errTime = Xh - solW_K(:,idxK,end);

    % Total error
    errTotal = Xh - solRef(:,end);

    errorTime_L2(i) = sqrt(mean(errTime.^2));
    errorTotal_L2(i) = sqrt(mean(errTotal.^2));
end

%% Estimate total convergence order

pTotal = polyfit(log(Kvalues),log(errorTotal_L2'),1);

fprintf('\n');
fprintf('====================================================\n');
fprintf('TOTAL ERROR\n');
fprintf('====================================================\n');
fprintf('Estimated total order with respect to K = %.4f\n', ...
        -pTotal(1));

%% Plot 3: total error

figure;

loglog(Kvalues,errorWZ_L2,'o-', ...
    'DisplayName','Wong-Zakai error ||X^K-X||');

hold on;

loglog(Kvalues,errorTime_L2,'s-', ...
    'DisplayName','time error ||X_h^K-X^K||');

loglog(Kvalues,errorTotal_L2,'^-', ...
    'DisplayName','total error ||X_h^K-X||');

grid on;

xlabel('number of Fourier modes K');
ylabel('L^2 error at T');

title('Error decomposition');

legend('Location','best');


%% ================================================================
% RK4 solver for the Wong-Zakai ODE
% ================================================================

function Xend = solve_WZ_RK4(X0, T, h, lambda, beta, ...
                               xiEndpoint, xiBridge, K)

    % Make T an exact endpoint
    N = ceil(T/h);
    h = T/N;

    X = X0;

    for n = 1:N

        t = (n-1)*h;

        k1 = rhsWZ( ...
            t,X,lambda,beta,xiEndpoint,xiBridge,K,T);

        k2 = rhsWZ( ...
            t+h/2,X+h*k1/2,...
            lambda,beta,xiEndpoint,xiBridge,K,T);

        k3 = rhsWZ( ...
            t+h/2,X+h*k2/2,...
            lambda,beta,xiEndpoint,xiBridge,K,T);

        k4 = rhsWZ( ...
            t+h,X+h*k3,...
            lambda,beta,xiEndpoint,xiBridge,K,T);

        X = X + h*(k1+2*k2+2*k3+k4)/6;

    end

    Xend = X;

end


%% ================================================================
% Right-hand side of Wong-Zakai ODE
% ================================================================

function f = rhsWZ(t,X,lambda,beta,...
                   xiEndpoint,xiBridge,K,T)

    % dWKdt = xiEndpoint/sqrt(T);
    % 
    % for k = 1:K
    %     dWKdt = dWKdt ...
    %         + sqrt(2/T)*xiBridge(k)*cos(k*pi*t/T);
    % end

    C_t = cos(pi/T * t * (1:K-1));
    dWKdt = xiEndpoint/sqrt(T) + sqrt(2/T)*C_t * xiBridge';
    
    f = lambda*X.*(1-X)-0.5*beta*X + beta*X*dWKdt;
end

% function dWKdt = dWZdt(t,T,xiEndpoint,xiBridge,K)
% 
%     dWKdt = xiEndpoint/sqrt(T);
% 
%     for k = 1:K
%         dWKdt = dWKdt ...
%             + sqrt(2/T)*xiBridge(k)*cos(k*pi*t/T);
%     end
% 
% end

function W = brownianFourierPath(t,xi,T)

    % Fourier representation
    %
    % W(t) =
    %
    %     xi_1 / sqrt(T) * t
    %
    %     + sqrt(2*T)/pi *
    %       sum_{k=2}^b xi_k/k * sin(k*pi*t/T)

    W = xi(1)/sqrt(T)*t;

    for k = 2:length(xi)

        W = W + ...
            sqrt(2*T)/pi * ...
            xi(k)/k .* sin(k*pi*t/T);

    end

end