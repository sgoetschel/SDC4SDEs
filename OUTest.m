%% ================================================================
%  Convergence study for OU process with Wong-Zakai approximation
%
%  SDE:
%
%      dX = lambda*(mu-X) dt + sigma dW (Ito and Stratonovich)
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
lambda = 2;
mu     = 1;
sigma  = 0.5;

nRealizations = 10000;
nStepsRef = 128;
tRef = linspace(0,T,nStepsRef+1);
dtRef = tRef(2)-tRef(1);

rng(12345);

%% ================================================================
% Fourier/Wong-Zakai representation of Brownian Bridge
% ================================================================
%
% We use
%
%   dW^K/dt =
%       xi_0/sqrt(T)
%       + sqrt(2/T) sum_{k=1}^K xi_k cos(k*pi*t/T)
%
% with xi_k ~ N(0,1).

Kvalues = [1 2 4 8 16];
Kmax = max(Kvalues);
nK = length(Kvalues);

% % Endpoint coefficient
% xiEndpoint = randn(nRealizations,1);
% 
% % Brownian bridge coefficients
% xiBridge = randn(nRealizations,Kmax);

% for pathwise/strong convergence, we need shared realizations
% construct sol from existing Brownian motion eta
eta = zeros(nRealizations, nStepsRef+1);
eta(:, 2:end) = cumsum(randn(nRealizations,nStepsRef),2).*sqrt(T)./sqrt(nStepsRef);
dW = eta(:,2:end)-eta(:,1:end-1);

xiEndpoint = eta(:,end);


% compute Wong-Zakai approximation *of the Brownian Bridge connecting 0 to xiEndpoint=W(T)*
% for k=2:Kmax    
%     phi = sqrt(2*T)/pi * sin(omega*tRef)/k;
%     an(:,k) = sum(phi .* eta,2);
% end
% Fourier coefficients
% n = (0:nStepsRef-1)';
k = 1:Kmax-1;                   % without the linear part
C = cos(pi/T * tRef(1:end-1)' * k); %cos(pi/nStepsRef * n * k);
xi = (sqrt(2/T) * C' * dW')';
% Reconstruction matrix
S = sin(pi/T * tRef' * k);
% Brownian-bridge coefficients
a = sqrt(2*T) ./ (pi*k) .* xi;


% check: re-construct Brownian motion from WZ
% All modes
W_Kmax = (tRef/T).*xiEndpoint + (S*a')';


% pathwise solution using given Brownian motion
% for r=1:nRealizations
%     Xpath(r,1) = X0;
%     for kk = 2:nStepsRef+1
%         dt = tRef(kk) - tRef(kk-1);
%         dW = eta(r,kk)-eta(r,kk-1); 
% 
%         e = exp(-lambda*dt);
%         m = mu + (Xpath(r,kk-1) - mu) * e;
%         v = (sigma^2 / (2*lambda)) * (1 - exp(-2*lambda*dt));
% 
%         Xpath(r,kk) = m + sqrt(max(v,0)) * dW/sqrt(dt);
%     end
%     XpathT(r) = Xpath(r,end);
% end
Xpath = zeros(nRealizations, nStepsRef+1);
Xpath(:,1) = X0;

for kk = 2:nStepsRef+1
    dt = tRef(kk) - tRef(kk-1);
    dW = eta(:,kk) - eta(:,kk-1);              % nRealizations x 1
    e  = exp(-lambda*dt);
    m  = mu + (Xpath(:,kk-1) - mu) * e;         % nRealizations x 1
    v  = (sigma^2 / (2*lambda)) * (1 - exp(-2*lambda*dt));
    Xpath(:,kk) = m + sqrt(max(v,0)) * dW/sqrt(dt);
end

XpathT = Xpath(:,end);

meanRef = mean(XpathT);
varRef  = var(XpathT);

mean_theory = mu + (X0 - mu)*exp(-lambda*T);
var_theory  = sigma^2/(2*lambda) * (1 - exp(-2*lambda*T));
meanKL = zeros(nK);
varKL = zeros(nK);
 
fprintf('\nReference statistics:\n');
fprintf('  E[X(T)]     = %.12e\t theory: %.12e\n',meanRef, mean_theory);
fprintf('  Var[X(T)]   = %.12e\t theory: %.12e\n',varRef, var_theory);

fprintf('\nSampling error: mean = %.3e\t var = %.3e\n', ...
    abs(mean_theory-meanRef), abs(var_theory-varRef));

% Check approximation of BB sol of OU wrt. modes
XK = zeros(nRealizations,nK);
XKfull = zeros(nRealizations,nK,nStepsRef+1);

deterministicPart = ...
    mu + (X0-mu)*exp(-lambda*tRef);

for i = 1:nK
    K = Kvalues(i);

    fprintf('\n----------------------------------------\n');
    fprintf('Number of Fourier modes: %d\n',K);

    W_K = (tRef/T).*xiEndpoint + (S(:,1:K-1)*a(:,1:K-1)')'; 

    % plot comparing one realization of Brownian motion and its BB-Fourier
    % approximation
    % figure;
    % plot(tRef, eta(1051,:), tRef, W_K(1051,:)) 

    % build OU solution for the current bridge
    % linear part first
    X_M = deterministicPart ...
        + sigma*(tRef/T).*xiEndpoint.*(1-exp(-lambda*tRef))/lambda;

    for kk = 1:K-1
        omega = kk*pi/T;
        I = ( ...
            lambda*cos(omega*tRef) ...
            + omega*sin(omega*tRef) ...
            - lambda*exp(-lambda*tRef) ...
            ) ./ (lambda^2 + omega^2);


        X_M = X_M ...
            + sigma*sqrt(2/T)*xi(:,kk)*I;
    end

    % plot paths corresponding to Brownian motion solution and BB-Fourier
    % figure;
    % plot(tRef, Xpath(1051,:), tRef, X_M(1051,:)) 


    errFull = X_M -Xpath;
    errorFullWZ_L2(i) = sqrt(mean(dtRef*sum(errFull.^2,2)));

    %XKfull(:,i,:) = X_M;
    XK(:,i) = X_M(:,end);
    meanBB(i) = mean(XK(:,i));
    varBB(i)  = var(XK(:,i));

    fprintf('KL mean %.12e\t KL var %.12e\n', meanBB(i), varBB(i));    
    fprintf('Weak mean error to theory = %.3e\n',abs(meanBB(i)-mean_theory));
    fprintf('Weak var  error to theory = %.3e\n',abs(varBB(i)-var_theory));
end

clearvars Xpath

% check variance error, can compute tail (variance due to missing Fourier
% modes) directly
omega = (1:Kmax)*pi/T;
IkT = lambda*((-1).^(1:Kmax) - exp(-lambda*T)) ...
      ./ (lambda^2 + omega.^2);

tailVar = zeros(nK-1,1);
for j = 1:nK-1
    tailVar(j) = (2*sigma^2/T) * sum(IkT(Kvalues(j)+1:end).^2);
end
% compute order for consecutive doubling of K, should be order 3
% tailVarRate = log2(tailVar(1:end-1)./tailVar(2:end));

% %% ================================================================
% % Coefficients for the OU stochastic integral
% % ================================================================
% %
% % The exact OU solution is
% %
% %   X(T) = mu + (X0-mu) exp(-lambda*T)
% %          + sigma * I,
% %
% % where
% %
% %   I = int_0^T exp(-lambda*(T-s)) dW(s).
% %
% % For the truncated Fourier representation,
% %
% %   I^K = sum_{k=1}^K a_k xi_k.
% %
% % These coefficients are known analytically.
% 
% a0 = (1-exp(-lambda*T))/(lambda*sqrt(T));
% 
% a = zeros(Kmax,1);
% 
% for k = 1:Kmax
% 
%     omega = k*pi/T;
% 
%     a(k) = sqrt(2/T) * ...
%         lambda*((-1)^k-exp(-lambda*T)) / ...
%         (lambda^2 + omega^2);
% 
% end
% 
% %% Variance of the exact OU stochastic integral
% 
% varI = (1-exp(-2*lambda*T))/(2*lambda);
% 
% %% ================================================================
% % Construct X^K(T) and exact X(T)
% % ================================================================
% 
% XK = zeros(nRealizations,nK);
% X  = zeros(nRealizations,nK);
% Xpath = zeros(nRealizations,nStepsRef+1);
% XpathT = zeros(nRealizations,1);
% 
% % Independent Gaussian variable representing the unresolved tail
% zTail = randn(nRealizations,1);
% 
% deterministicPart = ...
%     mu + (X0-mu)*exp(-lambda*T);
% 
% 
% 
% for i = 1:nK
% 
%     K = Kvalues(i);
% 
%     fprintf('\n----------------------------------------\n');
%     fprintf('Number of Fourier modes: %d\n',K);
% 
%     % Truncated stochastic integral
%     %IK = xi(:,1:K) * a(1:K);
%     IK = a0*xiEndpoint + xiBridge(:,1:K)*a(1:K);
% 
%     % ------------------------------------------------------------
%     % Exact Wong-Zakai solution
%     % ------------------------------------------------------------
% 
%     %XK(:,i) = deterministicPart + sigma*IK;
%     XK(:,i) = deterministicPart + sigma*IK;
% 
%     % ------------------------------------------------------------
%     % Exact OU solution
%     %
%     % I = I^K + I_tail
%     %
%     % I_tail is independent Gaussian with the variance of
%     % the unresolved Fourier modes.
%     % ------------------------------------------------------------
% 
%     varTail = varI - a0^2 - sum(a(1:K).^2);
% 
%     % Protect against tiny negative roundoff errors
%     varTail = max(varTail,0);
% 
%     Iexact = IK + sqrt(varTail)*zTail;
% 
%     X(:,i) = deterministicPart + sigma*Iexact;
% 
%     meanKL(i) = mean(XK(:,i));
%     varKL(i)  = var(XK(:,i));
% 
%     meanExKL(i) = mean(X(:,i));
%     varExKL(i)  = var(X(:,i));
% 
%     fprintf('KL mean %.12e\t KL var %.12e\n', meanKL(i), varKL(i));
%     %fprintf('ex mean %.12e\t ex var %.12e\n\n', meanExKL(i), varExKL(i));
% 
%     fprintf('Weak mean error  = %.3e\n',abs(meanKL(i)-mean_theory));
%     fprintf('Weak var error   = %.3e\n',abs(varKL(i)-var_theory));
%     fprintf('Weak mean error ex = %.3e\n',abs(meanExKL(i)-meanKL(i)));
%     fprintf('Weak var error  ex = %.3e\n',abs(varExKL(i)-varKL(i)));
% end

%% ================================================================
% 1. RK4 TIME-STEP CONVERGENCE
% ================================================================
%
% Fix K and vary h.
%
% The reference is X^K, which is known analytically.
%
% Therefore this measures ONLY:
%
%       X_h^K - X^K

Kfixed = 16;

hValues = T/(10*Kfixed) ./ 2.^(0:6);

nH = length(hValues);

errorTimeL1 = zeros(nH,1);
errorTimeL2 = zeros(nH,1);

idxK = find(Kvalues == Kfixed);

XKref = XK(:,idxK);

for ih = 1:nH

    h = hValues(ih);

    Xh = zeros(nRealizations,1);

    for r = 1:nRealizations

        Xh(r) = solveOU_WZ_RK4( ...
            X0,T,h,lambda,mu,sigma,...
            xiEndpoint(r), xi(r,1:Kfixed-1),Kfixed);
    end

    err = Xh - XKref;

    errorTimeL1(ih) = mean(abs(err));

    errorTimeL2(ih) = sqrt(mean(err.^2));

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

loglog(hValues,errorTimeL1,'o-','DisplayName','L^1 error');
hold on;
loglog(hValues,errorTimeL2,'s-','DisplayName','L^2 error');

grid on;

xlabel('time step h');
ylabel('error at T');

title(sprintf('RK4 time-step convergence, K = %d',Kfixed));

legend('Location','best');

%% ================================================================
% 2. WONG-ZAKAI CONVERGENCE
% ================================================================
%
% Compare the exact Wong-Zakai solution X^K with the exact
% Brownian OU solution X.
%
% Therefore:
%
%       error = X^K - X
%
% contains NO time discretization error.

errorWZ_L1 = zeros(nK,1);
errorWZ_L2 = zeros(nK,1);
errMean = zeros(nK,1);
errVar = zeros(nK,1);
errorVarOfError = zeros(nK,1);

for i = 1:nK

    err = XK(:,i)-XpathT;
    errorWZ_L1(i) = mean(abs(err));
    errorWZ_L2(i) = sqrt(mean(err.^2));

    errMean(i) = abs(meanBB(i)-mean_theory);
    errVar(i) = abs(varBB(i)-var_theory);
    errorVarOfError(i) = var(err,1);

    % % move this to above to avoid storing the full solution
    % errFull = squeeze(XKfull(:,i,:)) -Xpath;
    % errorFullWZ_L2(i) = sqrt(mean(dtRef*sum(errFull.^2,2)));
end

% Fit slopes with respect to K.
%
% If error ~ K^(-q), then
%
%   log(error) = -q log(K) + C.

pWZ_L1 = polyfit(log(Kvalues),log(errorWZ_L1'),1);
pWZ_L2 = polyfit(log(Kvalues),log(errorWZ_L2'),1);
pWZ_fullL2 = polyfit(log(Kvalues),log(errorFullWZ_L2'),1);

pWZ_weakMean = polyfit(log(Kvalues),log(errMean'),1);
pWZ_weakVar  = polyfit(log(Kvalues),log(errVar'),1);
pWZ_varOfErr  = polyfit(log(Kvalues),log(errorVarOfError'),1);
pWZ_tailVar = polyfit(log(Kvalues(2:end)),log(tailVar'),1);



fprintf('\n');
fprintf('====================================================\n');
fprintf('WONG-ZAKAI CONVERGENCE\n');
fprintf('====================================================\n');
fprintf('L1 order with respect to K = %.4f\n',-pWZ_L1(1));  % should be 3/2 at final time
fprintf('L2 order with respect to K = %.4f\n',-pWZ_L2(1));  % should be 3/2 at final time
fprintf('L2 order (all times) with respect to K = %.4f\n',-pWZ_fullL2(1)); % should be 1/2
fprintf('Weak (mean) order with respect to K = %.4f\n',-pWZ_weakMean(1));  % error in mean should be zero, so we only see sample error here
fprintf('Weak (var)  order with respect to K = %.4f\n',-pWZ_weakVar(1));   % this contains MCMC sampling error, probably different MCMC error from the next line. Should theoretically be the same, so order 3?
fprintf('Var of err  order with respect to K = %.4f\n',-pWZ_varOfErr(1));  % this should be order 3 as the tail variance
fprintf('Tail variance order with respect to K = %.4f\n',-pWZ_tailVar(1));  % this should be 3, it is already the expectation over the truncated Fourier coefficients


%% Plot 2: Wong-Zakai convergence

figure;

loglog(Kvalues,errorWZ_L1,'o-','DisplayName','Strong L^1 error');
hold on;
loglog(Kvalues,errorWZ_L2,'s-','DisplayName','Strong L^2 error');
loglog(Kvalues,errMean,'s-','DisplayName','Error in mean');
loglog(Kvalues,errVar,'s-','DisplayName','Error in variance');
loglog(Kvalues,errorVarOfError,'s-','DisplayName','Variance of error');
loglog(Kvalues,errorFullWZ_L2,'s-','DisplayName','L^2(0,T) error');


loglog(Kvalues, abs(mean_theory-meanRef)*ones(size(Kvalues)), 'DisplayName', 'Sampling error mean');
loglog(Kvalues, abs(var_theory-varRef)*ones(size(Kvalues)), 'DisplayName', 'Sampling error var')

grid on;

xlabel('number of Fourier modes K');
ylabel('error at T');

title('Wong-Zakai convergence');

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

        Xh(r) = solveOU_WZ_RK4( ...
            X0,T,h,lambda,mu,sigma,...
            xiEndpoint(r), xi(r,1:K-1),K);
       
    end

    % Time discretization component
    errTime = Xh - XK(:,i);

    % Total error
    errTotal = Xh - XpathT;

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

function Xend = solveOU_WZ_RK4(X0, T, h, lambda, mu, sigma, ...
                               xiEndpoint, xiBridge, K)

    % Make T an exact endpoint
    N = ceil(T/h);
    h = T/N;

    X = X0;

    for n = 1:N

        t = (n-1)*h;

        k1 = rhsWZ( ...
            t,X,lambda,mu,sigma,xiEndpoint,xiBridge,K,T);

        k2 = rhsWZ( ...
            t+h/2,X+h*k1/2,...
            lambda,mu,sigma,xiEndpoint,xiBridge,K,T);

        k3 = rhsWZ( ...
            t+h/2,X+h*k2/2,...
            lambda,mu,sigma,xiEndpoint,xiBridge,K,T);

        k4 = rhsWZ( ...
            t+h,X+h*k3,...
            lambda,mu,sigma,xiEndpoint,xiBridge,K,T);

        X = X + h*(k1+2*k2+2*k3+k4)/6;

    end

    Xend = X;

end

% function Xend = solveOU_WZ_RK4(X0, T, h, lambda, mu, sigma, ...
%                                xiEndpoint, xiBridge, K)
% 
%     % Make T an exact endpoint
%     N = ceil(T/h);
%     h = T/N;
% 
%     X = X0;
% 
%     for n = 1:N
% 
%         t = (n-1)*h;
% 
%         k1 = rhsWZ( ...
%             t,X,lambda,mu,sigma,xiEndpoint,xiBridge,K,T);
% 
%         k2 = rhsWZ( ...
%             t+h/2,X+h*k1/2,...
%             lambda,mu,sigma,xiEndpoint,xiBridge,K,T);
% 
%         k3 = rhsWZ( ...
%             t+h/2,X+h*k2/2,...
%             lambda,mu,sigma,xiEndpoint,xiBridge,K,T);
% 
%         k4 = rhsWZ( ...
%             t+h,X+h*k3,...
%             lambda,mu,sigma,xiEndpoint,xiBridge,K,T);
% 
%         X = X + h*(k1+2*k2+2*k3+k4)/6;
% 
%     end
% 
%     Xend = X;
% 
% end


%% ================================================================
% Right-hand side of Wong-Zakai ODE
% ================================================================

function f = rhsWZ(t,X,lambda,mu,sigma,...
                   xiEndpoint,xiBridge,K,T)
    
    C_t = cos(pi/T * t * (1:K-1));
    dWKdt = xiEndpoint/sqrt(T) + sqrt(2/T)*C_t * xiBridge';

    % dWKdt = xiEndpoint/T;
    % for k = 1:K-1
    %     dWKdt = dWKdt ...
    %          + sqrt(2/T)*xiBridge(k)*cos(k*pi*t/T); 
    % end

    f = lambda*(mu-X) + sigma*dWKdt;

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

