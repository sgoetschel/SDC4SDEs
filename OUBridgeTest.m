%% ================================================================
%  Convergence study for a nonlinear SDE with Wong-Zakai approximation
%
%  SDE:
%
%     dX = lambda*(mu-X) dt + sigma dW
%
%  Use KL-representation as exact for some number of modes,
%  for Delta T convergence, use same, or higher number of modes to
%  have exact representation of the noise
% ================================================================

clearvars;
close all;
clc;

%% Parameters
T      = 1;
X0     = 0.5;
lambda = 2.0;
mu     = 1.0;
sigma  = 0.5;

nStepsRef     = 128*1024;                                  % max number of time steps at which Brownian motion is sampled
nSteps        = [4, 8, 16, 32, 64, 128];      % time steps for convergence study
nStepsizes    = length(nSteps);
nRealizations = 100;

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

KRef = 64;   % "true" representation
KTest = 128;  % equal/larger than KRef allows to get correct Bridge, smaller is yet another approximation


%% ================================================================
% Analytical solution
% ================================================================

solRef = zeros(nRealizations, nStepsRef+1);
solRefnonsmooth = zeros(nRealizations, nStepsRef+1);

% sample Brownian motion at finest resolution in time
eta = zeros(nRealizations, nStepsRef+1);
eta(:, 2:end) = cumsum(randn(nRealizations,nStepsRef),2).*sqrt(T)./sqrt(nStepsRef);

dW = eta(:,2:end)-eta(:,1:end-1);  % increments of Brownian motion
xiEndpoint = eta(:,end);           % end value of Brownian motion after nStepsRef time steps


% compute reference Wong-Zakai approximation *of the Brownian Bridge connecting 0 to xiEndpoint=W(T)*
k = 1:KRef-1;                        % without the linear part
C = cos(pi/T * tRef(1:end-1)' * k);  %cos(pi/nStepsRef * n * k);
xi = (sqrt(2/T) * C' * dW')';
% Reconstruction matrix
S = sin(pi/T * tRef' * k);
% Brownian-bridge coefficients
a = sqrt(2*T) ./ (pi*k) .* xi;

W_KRef = zeros(nRealizations, nStepsRef+1);  % "true" driving Brownian motion
dW_KRef = zeros(nRealizations, nStepsRef); % ... and the "true" increments

%% ================================================================
%  Evaluate Wong-Zakai approximation error for the full solution
%  No use of Brownian bridges here!
%  Reference nonsmooth: with sampled Brownian motion on fine time grid
%  Ref Approximation: use KL-expansion instead of sampled Brownian motion

for  l=1:nRealizations
    W_KRef(l,:) = (tRef/T).*xiEndpoint(l) + (S(:,1:KRef-1)*a(l,1:KRef-1)')'; 
    dW_KRef(l,:) = W_KRef(l,2:end)-W_KRef(l,1:end-1);
    solRef(l,:) = exactOU([lambda, mu], sigma, 0, T, dtRef, X0, W_KRef(l,:));
    solRefnonsmooth(l,:) = exactOU([lambda, mu], sigma, 0, T, dtRef, X0, eta(l,:));
end

err_smooth_W_approx = solRefnonsmooth(:,end)-solRef(:,end);
errorSmoothL1 = mean(abs(err_smooth_W_approx));
errorSmoothL1weak = abs(mean(solRefnonsmooth(:,end))-mean(solRef(:,end)));
fprintf('W vs smooth W ref:\n\t strong error = %d\n\t weak error = %d\n',errorSmoothL1,errorSmoothL1weak);

meanRef = mean(solRef(:,end));
varRef  = var(solRef(:,end));
meanRefnonsmooth = mean(solRefnonsmooth(:,end));
varRefnonsmooth  = var(solRefnonsmooth(:,end));
trueMean = mu + (X0-mu)*exp(-lambda*T);
trueVar = sigma^2/(2*lambda)*(1-exp(-2*lambda*T));
fprintf('sample error vs true mean/var for nonsmooth W:\n\t mean = %d\n\t var = %d\n', abs(trueMean-meanRefnonsmooth), abs(trueVar-varRefnonsmooth));
fprintf('----------------------------\n');

%% ================================================================
%  Exact bridge convergence study
% ================================================================

K = KTest-1;   % without the linear part
combinedSol = zeros(nRealizations,nStepsRef+1);

% SDC setup
col_points = 5;
max_iter = 20;
nodes = 'lobatto';
tol = 1e-10;
strInit = 'constInit';
nComponents = 1;

[rhs, stochRhs, J, RhsIto, exact] = problemOU([lambda, mu], sigma);

allBridgeErrors = zeros(nStepsizes, 3);
allBridgenonsmoothErrors  = zeros(nStepsizes, 3);
allSDCErrors = zeros(nStepsizes, 3);
allSDCnonsmoothErrors = zeros(nStepsizes, 3);

rhsEvals = zeros(nStepsizes,1);
data =zeros(nStepsizes,14);        % data to store for plotting


for iN = 1:nStepsizes
    N = nSteps(iN);
    h = T/N;

    solAtTs = zeros(nRealizations,N+1);
    solSDCAtTs  = zeros(nRealizations,N+1);


    fprintf('N = %d\n',N);
    fprintf('----------------------------\n');
        
    % Number of fine reference increments per coarse interval
    nFinePerCoarse = nStepsRef/N;

    xi_l = zeros(nRealizations,K,N);
    
    for r = 1:nRealizations
        X0loc = X0;
        solWK_loc = zeros(N,nFinePerCoarse+1);

        solAtTs(r, 1) = X0loc;
        
        for j = 1:N
            idx = (j-1)*nFinePerCoarse + ...
                  (1:nFinePerCoarse);
            dWloc = dW_KRef(r, idx);
            DeltaW = sum(dWloc);
            xi = localKLcoefficients(dWloc,h,K);
            xi_l(r,:,j) = sqrt(2/h)*xi;

            t_subinterval_start = (j-1)*h;
            t_subinterval_end = j*h;
            t_subinterval_full = [tRef(idx) t_subinterval_end];

            Smat = sin(pi * ((t_subinterval_full'-t_subinterval_start)/h * (1:K)));

            % build local bridge
            etaLeft = W_KRef(r,(j-1)*nFinePerCoarse+1);
            etaRight = W_KRef(r,j*nFinePerCoarse+1);
            WK_loc = etaLeft + (t_subinterval_full-t_subinterval_start)/h.*(etaRight-etaLeft) + (Smat(:,1:K)*(sqrt(2*h)./ (pi*(1:K)).*xi')')'; 
            
            solWK_loc(j,:) = exactOU([lambda, mu], sigma, t_subinterval_start, t_subinterval_end, dtRef, X0loc, WK_loc);
            X0loc = solWK_loc(j,end);
            combinedSol(r,idx) = solWK_loc(j,1:end-1);
            solAtTs(r,j+1)=X0loc;
        end
        combinedSol(r,end) = solWK_loc(N,end);

        % run SDC
        eta0 = (W_KRef(r,1:nFinePerCoarse:end)- [0, W_KRef(r,1:nFinePerCoarse:end-1)]) ./ sqrt(h);
        parameters = [col_points, N, max_iter, 0, T];
        [Smat, tt, ~] = computeSpecMat(0, T, h, N, col_points, nodes);
        quadMatK_c = quadMatKFun(tt, h, nodes, K);
        [solSDCfull, countRhsEvaluationsThisRun, correctionNorms, solSDCAtTimesteps] = SDC_SDE_BB_NEW(...
            parameters, X0, Smat, quadMatK_c, tt, h, nodes, sigma, rhs, stochRhs, ...
            eta0, DeltaW, RhsIto, nComponents, xi_l(r,:,:), strInit, K+1, tol);
        solSDCAtTs(r,:) = solSDCAtTimesteps;
        rhsEvals(iN) = rhsEvals(iN)+countRhsEvaluationsThisRun;
    end

    err_SBB_to_Ref_at_T = (combinedSol(:, end) - solRef(:,end)); 
    errorSBBL1 = mean(abs(err_SBB_to_Ref_at_T));
    errorSBBL2 = sqrt(mean(err_SBB_to_Ref_at_T.^2));
    errWeakL1 = abs(mean(combinedSol(:,end))-mean(solRef(:,end)));
    allBridgeErrors(iN, :) = [errorSBBL1; errorSBBL2; errWeakL1];

    err_SBB_to_nonsmoothRef_at_T = (combinedSol(:, end) - solRefnonsmooth(:,end)); 
    errorSBBL1 = mean(abs(err_SBB_to_nonsmoothRef_at_T));
    errorSBBL2 = sqrt(mean(err_SBB_to_nonsmoothRef_at_T.^2));
    errWeakL1 = abs(mean(combinedSol(:,end))-mean(solRefnonsmooth(:,end)));
    allBridgenonsmoothErrors(iN, :) = [errorSBBL1; errorSBBL2; errWeakL1];

    err_SDC_to_Ref_at_T = (solSDCAtTs(:, end) - solRef(:,end)); 
    errorSDCL1 = mean(abs(err_SDC_to_Ref_at_T));
    errorSDCL2 = sqrt(mean(err_SDC_to_Ref_at_T.^2));
    errWeakL1 = abs(mean(solSDCAtTs(:,end))-mean(solRef(:,end)));
    allSDCErrors(iN, :) = [errorSDCL1; errorSDCL2; errWeakL1];

    err_SDC_to_nonsmoothRef_at_T = (solSDCAtTs(:, end) - solRefnonsmooth(:,end)); 
    errorSDCL1 = mean(abs(err_SDC_to_nonsmoothRef_at_T));
    errorSDCL2 = sqrt(mean(err_SDC_to_nonsmoothRef_at_T.^2));
    errWeakL1 = abs(mean(solSDCAtTs(:,end))-mean(solRefnonsmooth(:,end)));
    allSDCnonsmoothErrors(iN, :) = [errorSDCL1; errorSDCL2; errWeakL1];

    % store diagnostics data
    data(iN, :) = [h, col_points, KRef, KTest, nRealizations, rhsEvals(iN), ...
                                  allBridgeErrors(iN, 1), allBridgeErrors(iN, 3), ...
                                  allBridgenonsmoothErrors(iN, 1), allBridgenonsmoothErrors(iN, 3), ...
                                  allSDCErrors(iN, 1), allSDCErrors(iN, 3), ...
                                  allSDCnonsmoothErrors(iN, 1), allSDCnonsmoothErrors(iN, 3)];
end

% write data file
filename = sprintf('ou_data_%s.dat', string(datetime('now','TimeZone','local','Format','d-MMM-y-HH:mm:ss')));
fileID = fopen(filename,'w');
fprintf(fileID, ['stepSize \t colpoints \t KRef \t KTest \t nRealizations \t rhsEvalsSum \t ' ...
    ' SBB-to-Ref-strong \t SBB-to-Ref-weak \t ' ...
    ' SBB-to-nonsm-Ref-strong \t SBB-to-nonsm-Ref-weak \t' ...
    ' SBB-SDC-to-Ref-strong \t SBB-SDC-to-Ref-weak \t' ...
    ' SBB-SDC-to-nonsm-Ref-strong \t SBB-SDC-to-nonsm-Ref-weak \n']);
fprintf(fileID, ['%f \t %d \t %d \t %d \t %d \t %d \t' ...
    ' %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n'], data');
fclose(fileID);

% Estimate slopes
hValues = T./nSteps;
% Reference slopes dt^1, dt^2
C1 = allSDCErrors(1, 1) * hValues(1)^(-1);
C2 = allSDCErrors(1, 1) * hValues(1)^(-2);
C3 = allSDCErrors(1, 1) * hValues(1)^(-3);
C4 = allSDCErrors(1, 1) * hValues(1)^(-4);

fprintf('----------------------------\n');
fprintf('Strong convergence orders:\n');
fprintf('----------------------------\n');

%pSBBTimeL1 = polyfit(log(hValues),log(allBridgeErrors(:,1)'),1);
%pSBBTimeL2 = polyfit(log(hValues),log(allBridgeErrors(:,2)'),1);
%fprintf('SBB vs smooth W ref, temporal L1 order = %.4f\n',pSBBTimeL1(1));
%fprintf('SBB vs smooth W ref, temporal L2 order = %.4f\n',pSBBTimeL2(1));  % should be constant

pSDCTimeL1 = polyfit(log(hValues),log(allSDCErrors(:,1)'),1);
pSDCTimeL2 = polyfit(log(hValues),log(allSDCErrors(:,2)'),1);
fprintf('SDC-SBB vs smooth W ref, temporal L1 order = %.4f\n',pSDCTimeL1(1));
fprintf('SDC-SBB vs smooth W ref, temporal L2 order = %.4f\n',pSDCTimeL2(1));

pSDCTimeL1 = polyfit(log(hValues),log(allSDCnonsmoothErrors(:,1)'),1);
pSDCTimeL2 = polyfit(log(hValues),log(allSDCnonsmoothErrors(:,2)'),1);
fprintf('SDC-SBB vs nonsmooth W ref, temporal L1 order = %.4f\n',pSDCTimeL1(1));
fprintf('SDC-SBB vs nonsmooth W ref, temporal L2 order = %.4f\n',pSDCTimeL2(1));

fprintf('----------------------------\n');
fprintf('Weak convergence orders:\n');
fprintf('----------------------------\n');
%pSBBTimeL1 = polyfit(log(hValues),log(allBridgeErrors(:,3)'),1);
%fprintf('SBB vs smooth W ref, temporal L1 order = %.4f\n',pSBBTimeL1(1));  % should be constant

pSDCTimeL1 = polyfit(log(hValues),log(allSDCErrors(:,3)'),1);
fprintf('SDC-SBB vs smooth W ref, temporal L1 order = %.4f\n',pSDCTimeL1(1));

pSDCTimeL1 = polyfit(log(hValues),log(allSDCnonsmoothErrors(:,3)'),1);
fprintf('SDC-SBB vs nonsmooth W ref, temporal L1 order = %.4f\n',pSDCTimeL1(1));


set(0,'DefaultLineLineWidth',2)
%% Plot 1: strong convergence
figure;
hold on;
loglog(hValues,allBridgeErrors(:,1),'o-','DisplayName','SBB L^1 error');
loglog(hValues,allBridgenonsmoothErrors(:,1),'o-','DisplayName','SBB to nonsmooth L^1 error');
loglog(hValues,allSDCErrors(:,1),'o-','DisplayName','SDC-SBB L^1 error');
loglog(hValues,allSDCnonsmoothErrors(:,1),'o-','DisplayName','SDC-SBB to nonsmooth L^1 error');
loglog(hValues,C1*hValues.^(1),'--', 'DisplayName', 'slope 1');
loglog(hValues,C2*hValues.^(2),'--', 'DisplayName', 'slope 2');
loglog(hValues,C3*hValues.^(3),'--', 'DisplayName', 'slope 3');
loglog(hValues,C4*hValues.^(4),'--', 'DisplayName', 'slope 4');
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

grid on;
xlabel('time step h');
ylabel('error at T');
title(sprintf('strong convergence, K = %d',KTest));
legend('Location','best');

%% Plot 2: weak convergence
C1 = allSDCErrors(1, 3) * hValues(1)^(-1);
C2 = allSDCErrors(1, 3) * hValues(1)^(-2);
C3 = allSDCErrors(1, 3) * hValues(1)^(-3);

figure;
hold on;
loglog(hValues,allBridgeErrors(:,3),'o-','DisplayName','SBB L^1 error');
loglog(hValues,allBridgenonsmoothErrors(:,3),'o-','DisplayName','SBB to nonsmooth L^1 error');
loglog(hValues,allSDCErrors(:,3),'o-','DisplayName','SDC-SBB L^1 error');
loglog(hValues,allSDCnonsmoothErrors(:,3),'o-','DisplayName','SDC-SBB to nonsmooth L^1 error');
loglog(hValues,C1*hValues.^(1),'--', 'DisplayName', 'slope 1');
loglog(hValues,C2*hValues.^(2),'--', 'DisplayName', 'slope 2');
loglog(hValues,C3*hValues.^(3),'--', 'DisplayName', 'slope 3');
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

grid on;
xlabel('time step h');
ylabel('error at T');
title(sprintf('weak convergence, K = %d',KTest));
legend('Location','best');

% 
% %% ================================================================
% %  Exact bridge convergence study
% % ================================================================
% %
% % For each coarse grid and M, compute:
% %
% %     X_{h,M}^{exact}(T)
% %
% % and compare with Xref.
% %
% % This measures the bridge/Wong-Zakai error without RK4 error.
% 
% fprintf('\n');
% fprintf('====================================================\n');
% fprintf('BRIDGE CONVERGENCE\n');
% fprintf('====================================================\n');
% 
% errorBridge = zeros(nStepsizes,nK);
% solW_K(nRealizations,i,:) = exactOU([lambda, mu], sigma, 0, T, dtRef, X0, W_K);
% 
% errorSBBL1 = zeros(nK,1);
% errorSBBL2 = zeros(nK,1);
% errorSDCL1 = zeros(nK,1);
% errorSDCL2 = zeros(nK,1);
% errorSDCRefL1 = zeros(nK,1);
% errorSDCRefL2 = zeros(nK,1);
% 
% 
% allErrors = zeros(nStepsizes, nK, nRealizations);
% allSDCErrors = zeros(nStepsizes, nK, nRealizations);
% allSDCRefErrors = zeros(nStepsizes, nK, nRealizations);
% 
% 
% for iN = 1:nStepsizes
%     combinedSol = zeros(nRealizations,nK,nStepsRef+1);
% 
% 
%     N = nSteps(iN);
%     h = T/N;
% 
%     solAtTs = zeros(nRealizations,N+1);
%     solSDCAtTs  = zeros(nRealizations,N+1);
% 
%     fprintf('N = %d\n',N);
%     fprintf('----------------------------\n');
% 
%     % Number of fine reference increments per coarse interval
%     nFinePerCoarse = nStepsRef/N;
% 
%     % Compute the solution for every K independently.
%     %
%     for iK = 1:nK
%         K = Kvalues(iK);
% 
%         xi_l = zeros(nRealizations,K,N);
% 
% 
%         for r = 1:nRealizations
%             X0loc = X0;
%             solWK_loc = zeros(N,nFinePerCoarse+1);
% 
%             solAtTs(r, 1) = X0loc;
% 
%             for j = 1:N
%                 idx = (j-1)*nFinePerCoarse + ...
%                       (1:nFinePerCoarse);
%                 dWloc = dW(r, idx);
%                 DeltaW = sum(dWloc);
%                 xi = localKLcoefficients(dWloc,h,K);
%                 xi_l(r,:,j) = sqrt(2/h)*xi;
%                 t_subinterval = tRef(idx); % without endpoint
%                 t_subinterval_end = tRef(idx(end)+1);
%                 t_subinterval_full = [t_subinterval t_subinterval_end];
%                 t_subinterval_start = t_subinterval(1);
%                 Smat = sin(pi * ((t_subinterval_full'-t_subinterval_start)/h * (1:K)));
% 
%                 % build local bridge
%                 etaLeft = eta(r,(j-1)*nFinePerCoarse+1);
%                 etaRight = eta(r,j*nFinePerCoarse+1);
%                 WK_loc = etaLeft + (t_subinterval_full-t_subinterval_start)/h.*(etaRight-etaLeft) + (Smat(:,1:K)*(sqrt(2*h)./ (pi*(1:K)).*xi')')'; 
%                 solWK_loc(j,:) = exactOU([lambda, mu], sigma, t_subinterval_start, t_subinterval_end, dtRef, X0loc, WK_loc);
%                 X0loc = solWK_loc(j,end);
%                 combinedSol(r,iK,idx) = solWK_loc(j,1:end-1);
%                 solAtTs(r,j+1)=X0loc;
%             end
%             combinedSol(r,iK,end) = solWK_loc(N,end);
% 
%             %eta0 = eta(r,1:nFinePerCoarse:end) ./ sqrt(h); % divide by sqrt(h) for stupid reasons (SDC divides by sqrt(h), in total need to divide by h)
%             eta0 = (eta(r,1:nFinePerCoarse:end)- [0, eta(r,1:nFinePerCoarse:end-1)]) ./ sqrt(h);
%             parameters = [col_points, N, max_iter, t_begin, t_end];
%             [S, tt, ~] = computeSpecMat(t_begin, t_end, h, N, col_points, nodes);
%             quadMatK_c = quadMatKFun(tt, h, nodes, K);
%             [solSDCfull, countRhsEvaluationsThisRun, correctionNorms, solSDCAtTimesteps] = SDC_SDE_BB_NEW(...
%                 parameters, X0, S, quadMatK_c, tt, h, nodes, sigma, rhs, stochRhs, ...
%                 eta0, DeltaW, RhsIto, nComponents, xi_l(r,:,:), strInit, K+1, tol);
%             solSDCAtTs(r,:) = solSDCAtTimesteps;
%         end
%         err_to_Ref = (combinedSol(:,iK, end) - solRef(:,end)); % solRef(:,end) - solSDCAtTs(:,end); %
%         errorSBBL1(iK) = mean(abs(err_to_Ref));
%         errorSBBL2(iK) = sqrt(mean(err_to_Ref.^2));
%         allErrors(iN, iK, :) = err_to_Ref;
% 
%         err_SDC_to_Ref = combinedSol(:,iK, end) - solSDCAtTs(:,end);
%         errorSDCL1(iK) = mean(abs(err_SDC_to_Ref));
%         errorSDCL2(iK) = sqrt(mean(err_SDC_to_Ref.^2));
%         allSDCErrors(iN, iK, :) = err_SDC_to_Ref;
% 
%         err_SDC_to_Ref = solRef(:,end) - solSDCAtTs(:,end);
%         errorSDCRefL1(iK) = mean(abs(err_SDC_to_Ref));
%         errorSDCRefL2(iK) = sqrt(mean(err_SDC_to_Ref.^2));
%         allSDCRefErrors(iN, iK, :) = err_SDC_to_Ref;
%         % ttt = linspace(0,T,N+1);
%         % figure;
%         % hold on;
%         % plot(ttt, solAtTs, 'DisplayName', 'exact');
%         % plot(ttt, solSDCAtTimesteps, 'DisplayName', 'SDC');
%         % title(sprintf('N = %d, K = %d',N,K));
%         % legend;
%     end
%     % Estimate slopes
%     pSBBTimeL1 = polyfit(log(Kvalues),log(errorSBBL1'),1);
%     pSBBTimeL2 = polyfit(log(Kvalues),log(errorSBBL2'),1);
%     fprintf('L1 order = %.4f\n',-pSBBTimeL1(1));
%     fprintf('L2 order = %.4f\n',-pSBBTimeL2(1));
% 
%     figure;
%     loglog(Kvalues,errorSBBL1,'o-','DisplayName','L^1 error to ref');
%     hold on;
%     loglog(Kvalues,errorSBBL2,'s-','DisplayName','L^2 error to ref');
% 
%     grid on;
% 
%     xlabel('K');
%     ylabel('error at T');
%     title(sprintf('Bridge convergence, time-steps = %d',N));
%     legend('Location','best');
% end
% 
% timeErrL1 = zeros(nK,nStepsizes);
% timeErrL2 = zeros(nK,nStepsizes);
% figure;
% hold on;
% for iK=1:nK
%     timeErrL1(iK,:) = mean(abs(allErrors(:,iK,:)),3);
%     timeErrL2(iK,:) = mean(allErrors(:,iK,:).^2,3);
%     %loglog(nSteps,squeeze(timeErrL2(iK,:)), 'DisplayName', sprintf('L2, K=%d',Kvalues(iK)));
%     loglog(T./nSteps,squeeze(timeErrL1(iK,:)), 'DisplayName', sprintf('SBB+Exact to Ref L1, K=%d',Kvalues(iK)));
% end
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% legend('Location', 'best');
% 
% timeErrSDCL1 = zeros(nK,nStepsizes);
% timeErrSDCL2 = zeros(nK,nStepsizes);
% figure;
% hold on;
% for iK=1:nK
%     timeErrSDCL1(iK,:) = mean(abs(allSDCErrors(:,iK,:)),3);
%     timeErrSDCL2(iK,:) = mean(allSDCErrors(:,iK,:).^2,3);
%     %loglog(nSteps,squeeze(timeErrL2(iK,:)), 'DisplayName', sprintf('L2, K=%d',Kvalues(iK)));
%     loglog(T./nSteps,squeeze(timeErrSDCL1(iK,:)), 'DisplayName', sprintf('SBB+SDC to SBB L1, K=%d',Kvalues(iK)),'LineWidth',2);
% 
%     % Reference slope dt^(1/2)
%     C = timeErrSDCL1(iK,1) * (T./nSteps(1))^(-1);
% 
%     loglog((T./nSteps),C*(T./nSteps).^(1),'--',...
%         'LineWidth',1.2);
% end
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% legend('Location', 'best');
% 
% figure;
% hold on;
% for iK=1:nK
%     timeErrSDCL1(iK,:) = mean(abs(allSDCRefErrors(:,iK,:)),3);
%     timeErrSDCL2(iK,:) = mean(allSDCRefErrors(:,iK,:).^2,3);
%     %loglog(nSteps,squeeze(timeErrL2(iK,:)), 'DisplayName', sprintf('L2, K=%d',Kvalues(iK)));
%     loglog(T./nSteps,squeeze(timeErrSDCL1(iK,:)), 'DisplayName', sprintf('SBB+SDC to Ref L1, K=%d',Kvalues(iK)),'LineWidth',2);
% 
%     % Reference slope dt^(1/2)
%     C = timeErrSDCL1(iK,1) * (T./nSteps(1))^(-1);
% 
%     loglog((T./nSteps),C*(T./nSteps).^(1),'--',...
%         'LineWidth',1.2);
% end
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% legend('Location', 'best');
% 
% %% ================================================================
% % RK4 TIME-STEP CONVERGENCE
% % ================================================================
% %
% % Fix K and vary h.
% %
% % The reference is X^K, which is known analytically.
% %
% % Therefore this measures ONLY:
% %
% %       X_h^K - X^K
% 
% Kfixed = 64;
% 
% hValues = T/(Kfixed) ./ 2.^(0:6);
% 
% nH = length(hValues);
% 
% errorTimeL1 = zeros(nH,1);
% errorTimeL2 = zeros(nH,1);
% 
% idxK = find(Kvalues == Kfixed);
% 
% for ih = 1:nH
% 
%     h = hValues(ih);
% 
%     Xh = zeros(nRealizations,1);
% 
%     for r = 1:nRealizations
% 
%         Xh(r) = solve_WZ_RK4( ...
%             X0,T,h,lambda,beta,...
%             xiEndpoint(r), xi(r,1:Kfixed-1),Kfixed);
%     end
% 
%     err_to_Ref = Xh - solRef(:,end);
%     err  = Xh - solW_K(:,idxK,end);
% 
%     errorTimeL1(ih) = mean(abs(err));
% 
%     errorTimeL2(ih) = sqrt(mean(err.^2));
% 
%     errorRefTimeL1(ih) = mean(abs(err_to_Ref));
%     errorRefTimeL2(ih) = sqrt(mean(err_to_Ref.^2));
% end
% 
% % Estimate slopes
% pSBBTimeL1 = polyfit(log(hValues),log(errorTimeL1'),1);
% pSBBTimeL2 = polyfit(log(hValues),log(errorTimeL2'),1);
% 
% fprintf('\n');
% fprintf('====================================================\n');
% fprintf('RK4 TIME-STEP CONVERGENCE\n');
% fprintf('====================================================\n');
% fprintf('K = %d\n',Kfixed);
% fprintf('L1 order = %.4f\n',pSBBTimeL1(1));
% fprintf('L2 order = %.4f\n',pSBBTimeL2(1));
% 
% %% Plot 1: time-step convergence
% 
% figure;
% 
% loglog(hValues,errorTimeL1,'o-','DisplayName','L^1 error to WZ');
% hold on;
% loglog(hValues,errorTimeL2,'s-','DisplayName','L^2 error to WZ');
% loglog(hValues,errorRefTimeL1,'s-','DisplayName','L^1 error to ref');
% loglog(hValues,errorRefTimeL2,'s-','DisplayName','L^2 error to ref');
% 
% grid on;
% 
% xlabel('time step h');
% ylabel('error at T');
% 
% title(sprintf('RK4 time-step convergence, K = %d',Kfixed));
% 
% legend('Location','best');
% 
% %% ================================================================
% % 3. TOTAL ERROR
% % ================================================================
% %
% % Now solve the Wong-Zakai ODE numerically:
% %
% %       X_h^K
% %
% % and compare directly against the exact OU solution:
% %
% %       X.
% %
% % Thus:
% %
% %       total error = X_h^K - X.
% %
% % We choose h proportional to 1/K and sufficiently small such
% % that the time discretization error is small compared with the
% % Wong-Zakai error.
% 
% C = 80;
% 
% errorTotal_L2 = zeros(nK,1);
% errorTime_L2 = zeros(nK,1);
% 
% for i = 1:nK
% 
%     K = Kvalues(i);
% 
%     % Resolve the highest Fourier frequency increasingly well
%     h = T/(C*K);
% 
%     Xh = zeros(nRealizations,1);
% 
%     for r = 1:nRealizations
% 
%         Xh(r) = solve_WZ_RK4( ...
%             X0,T,h,lambda,beta,...
%             xiEndpoint(r), xi(r,1:K-1),K);
% 
%     end
% 
%     % Time discretization component
%     errTime = Xh - solW_K(:,idxK,end);
% 
%     % Total error
%     errTotal = Xh - solRef(:,end);
% 
%     errorTime_L2(i) = sqrt(mean(errTime.^2));
%     errorTotal_L2(i) = sqrt(mean(errTotal.^2));
% end
% 
% %% Estimate total convergence order
% 
% pTotal = polyfit(log(Kvalues),log(errorTotal_L2'),1);
% 
% fprintf('\n');
% fprintf('====================================================\n');
% fprintf('TOTAL ERROR\n');
% fprintf('====================================================\n');
% fprintf('Estimated total order with respect to K = %.4f\n', ...
%         -pTotal(1));
% 
% %% Plot 3: total error
% 
% figure;
% 
% loglog(Kvalues,errorWZ_L2,'o-', ...
%     'DisplayName','Wong-Zakai error ||X^K-X||');
% 
% hold on;
% 
% loglog(Kvalues,errorTime_L2,'s-', ...
%     'DisplayName','time error ||X_h^K-X^K||');
% 
% loglog(Kvalues,errorTotal_L2,'^-', ...
%     'DisplayName','total error ||X_h^K-X||');
% 
% grid on;
% 
% xlabel('number of Fourier modes K');
% ylabel('L^2 error at T');
% 
% title('Error decomposition');
% 
% legend('Location','best');


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

function xi = localKLcoefficients(dW,h,M)
% Compute local Brownian-bridge KL coefficients from Brownian
% increments on an interval of length h.
%
% dW contains increments of a fixed fine Brownian path.
%
% xi_k =
% sqrt(2/h) sum_n cos(k*pi*tau_n/h) dW_n

    nSub = length(dW);

    n = (0:nSub-1)';

    k = 1:M;

    C = cos(pi/nSub * n*k);

    xi = sqrt(2/h) * (C' * dW');

end