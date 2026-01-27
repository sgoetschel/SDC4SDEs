% tol - SDC sweep tolerance for correction norm
% NNfinest - number of time steps in finest time discretization (for convergence order studies)
% step_size - time step size
% t_begin, t_end - time interval
% col_points - number of collocation points
% max_iter - maximum number of SDC sweeps
% realIter - number of samples (realizations)
% initialValue - initial value of SDE
% nodes - string specifying type of collocation nodes
% rhs - function handle for drift term (Stratonovic)
% stochRhs - function handle for diffusion term
% RhsIto - function handle for drift term in Ito form (for some initializations)
% beta - coefficients of independent Wiener processes
% order - dimension of SDE system (number of variables/equations)
% eta - random variable specifying Wiener process/bridge (linear part)
% xi - random variable specifying Wiener process/bridge (cosines)
% strInit - string specifying initialization method
% nBridgeTerms - number of terms in the KL-expansion

function [solsSDC, compTime] = sdeSolverSBBSDC(tol, NNfinest, step_size, t_begin, t_end, col_points, max_iter, realIter, initialValue, nodes, beta, rhs, stochRhs, RhsIto, order, eta, xi, strInit, nBridgeTerms)

intervals = round((t_end-t_begin)/step_size); % number of time intervals in discretization

%% setup SDC method
parameters = [col_points, intervals, max_iter, t_begin, t_end];
    
%computes spectral matrix (S), the micro time nodes (t) and macro time nodes (x)
[S, t, ~] = computeSpecMat(t_begin, t_end, step_size, intervals, col_points, nodes);
        
NNCurrent = intervals; %(t_end-t_begin) / step_size;
        
compTime = 0;
        
xi_l = zeros(order, nBridgeTerms-1, intervals);

solsSDC = zeros(order, realIter);

for l=1:realIter
    for k=1:order
        etaMat = eta{k}(l,:);%eta{k}(l,:); %always use l=1 for a tied down process? i.e. all going through the same nodes
        % only xi, i.e., the cosine parts in the bridge vary
        etaCurrSteps(k,:) = [0 etaMat(NNfinest/NNCurrent:NNfinest/NNCurrent:end)];
        eta0(k,:) = etaCurrSteps(k,1:end)- [0, etaCurrSteps(k,1:end-1)];
        if nBridgeTerms > 1
            xi_l(k,:,:) = xi{k}(l,:,NNfinest/NNCurrent:NNfinest/NNCurrent:end);
        end
     end
     eta0 = sqrt(1/step_size)*eta0; % eta0 now has variance 1 (before: variance deltaT)
            
     if(nBridgeTerms>1)
        %compute quadrature matrix for each expansion term
        [quadMatK_c] = quadMatKFun(t, step_size, nodes, nBridgeTerms-1);
     else
        quadMatK_c = 0;
     end
                                
     %dW approximation
     deltaW = etaCurrSteps(:,2:end)- etaCurrSteps(:,1:end-1);
            
     %calls SDC_SDE_BB_NEW.m to compute SDC approximation
     tic
     [solSDC, countRhsEvaluations] = SDC_SDE_BB_NEW(parameters, initialValue, S, quadMatK_c, t,step_size, nodes, beta, rhs, stochRhs, eta0, deltaW, RhsIto, order, xi_l, strInit, nBridgeTerms, tol);
     elapsedTime = toc;
     compTime = compTime +elapsedTime;
     solsSDC(:,l) =  solSDC(:,end);    
end %for loop realIter
    
end %function
