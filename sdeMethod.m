%% This function chooses between EM, Milstein and SDC method
function [errWeak, errStrong, errL2, errT, compTime, countRhsEvaluations] = sdeMethod(sde_solver, tol, NNfinest, step_size, t_begin, t_end, intervals, col_points, max_iter, realIter, initial, nodes, lambda, beta, rhs, stochRhs, RhsIto, J, order, eta, exact, d, xi, plot_Error, plot_Sol, strInit, m, solRefAll)
time = t_begin:step_size:t_end;
nTime = intervals + 1;
NNCurrent = (t_end-t_begin) / step_size;
downsampleFactor = NNfinest / NNCurrent;
countRhsEvaluations = 0;

%% choose SDC method
if strcmp(sde_solver,'SDC_BB')
    
    parameters = [col_points, intervals, max_iter, t_begin, t_end];
    
    %computes spectral matrix (S), the micro time nodes (t) and macro time nodes (x)
    [S, t, ~] = computeSpecMat(t_begin, t_end, step_size, intervals, col_points, nodes);

        sumSol=zeros(order, nTime);
        sumSolRef=zeros(order, nTime);
        diffSol = zeros(order, nTime);
        diff = 0;
        
        sumSolT2 = 0;
        
        compTime = 0;
        solsSDC = zeros(realIter, 1);
        etaCurrSteps = zeros(order, nTime);
        eta0 = zeros(order, nTime);
        xi_l = zeros(order, m-1, intervals);
        
        
        for l=1:realIter
            for k=1:order
                etaMat = eta{k}(l,:);%eta{k}(l,:); %always use l=1 for a tied down process? i.e. all going through the same nodes
                % only xi, i.e., the cosine parts in the bridge vary
                etaCurrSteps(k,:) = [0 etaMat(downsampleFactor:downsampleFactor:end)];
                eta0(k,:) = etaCurrSteps(k,1:end)- [0, etaCurrSteps(k,1:end-1)];
                xi_l(k,:,:) = xi{k}(l,:,downsampleFactor:downsampleFactor:end);
            end
            eta0 = sqrt(1/step_size)*eta0;
            % eta0 now has variance 1 (before: variance deltaT)

            %global thisRealization;
            thisRealization = l;
            
            if(m>1)
                %compute quadrature matrix for each expansion term
                [quadMatK_c] = quadMatKFun(t, step_size, nodes, m-1);
            else
                quadMatK_c = 0;
                %quadMatK_s = 0;
            end
                                
            %dW approximation
            deltaW = etaCurrSteps(:,2:end)- etaCurrSteps(:,1:end-1);
            
            %calls SDC_SDE_BB.m to compute SDC approximation
            tic
            [solSDC, countRhsEvaluationsThisRun] = SDC_SDE_BB(parameters, initial(), S, quadMatK_c, t,step_size, nodes, beta, rhs, stochRhs, eta0, deltaW, RhsIto, order, xi_l, strInit, m, tol);
            elapsedTime = toc;
            compTime = compTime +elapsedTime;
            countRhsEvaluations = countRhsEvaluations + countRhsEvaluationsThisRun;
            solSDC =  solSDC(:,1:col_points-1:end);
            
            %exact solution at macro time steps
%             if strcmp(d, 'TP3')
%                 steps = step_size;
%                 step_size = 1/NNfinest;
%                 etaCurr = [0 etaMat];
%                 solRef = exact(lambda, beta, t_begin, t_end, step_size, initial, etaCurr);
%                 solRef = solRef(1:NNfinest/NNCurrent:end);
%                 step_size = steps;
%                 
%             else
%                steps = step_size;
%                 %sol of the Stratonovich SDE
%                 etaCurr = [0 etaMat];
%                 solRef = exact(lambda, beta, t_begin, t_end, 1/NNfinest, initial, etaCurr);
%                 solRef = solRef(1:NNfinest/NNCurrent:end);
                
                %sol of the smooth Brownian Bridge ODE
                % on finest time grid
%                 steps = 1.0/NNfinest;
%                 for k=1:order
%                 etaFin(k,:) = [0 etaMat];
%                 etaFinest(k,:) = etaFin(k,1:end)- [0, etaFin(k,1:end-1)];
%                 end
%                 etaFinest = sqrt(1/steps)*etaFinest;
%                 k=1;
%                 bM(:,k) = brownianBridge(etaFinest(:,k+1), steps , steps, xiFinest);
%                 for k=2:NNfinest
%                    bM(:,k) =  brownianBridge(etaFinest(:,k+1), steps , steps, xiFinest)+bM(:,k-1);
%                 end
% 
%                   bbM(:,:) = [zeros(order,1) bM(:,:)];
%                   solRef = exactExp_BB(lambda, beta, t_begin, t_end, steps, initial, bbM);
%                   solRef = solRef(1:NNfinest/NNCurrent:end);
               
                               
%               %ref sol using matlab ode solver  
%                 y0 = initial;
%                 tsolAppro = [t_begin];
%                 solAppro = [initial];
%                 for k=1:intervals
%                     [tsol, solA] = ode45(@(s,x)SBB_ODE(s,x, eta0(:,k+1), steps, xi_l(:,:,k), lambda, beta),[t_begin+(k-1)*steps t_begin+k*steps], y0');
%                     y0 = solA(end,:);
%                     tsolAppro = [tsolAppro tsol(end)'];
%                     solAppro = [solAppro solA(end)'];
%                 end
%                 solRef = solAppro;

               % step_size = steps;
              if strcmp(d, 'OU') || strcmp(d, 'Mattingly')
                solRef= 0;
              else
               solRef = solRefAll(:,:,l);
              end
  %          end
            
            %adding values of each realization realIter
            %preperatory work for determining the expected solution
            sumSol = sumSol + solSDC;
            if strcmp(d, 'OU') || strcmp(d, 'Mattingly')
              sumSolT2 = sumSolT2 + solSDC(:,end).*solSDC(:,end);
            end
            sumSolRef = sumSolRef + solRef;
            
            %preparatory work for computation of different error definitions
            diffSol = diffSol + abs(solSDC-solRef);
            diff = diff + (solSDC(:,end)-solRef(:,end)).^2;
            solsSDC(l) = solSDC(1,end);
            
        end %for loop realIter
        
        %calculate expected solution for each time point t
        sol= sumSol/realIter;
        solExactRef = sumSolRef/realIter;
        
%         if strcmp(d, 'TP3') % compute reference solution new with sample paths
%           nstepsRef = 1000;
%           refdt = step_size/nstepsRef;
%           solsRef = zeros(realIter,1);
%           reft = refdt:refdt:step_size;
%           for l=1:realIter
%             refWt = cumsum(randn(1,nstepsRef)).*sqrt(refdt); 
%             %delta=refWt(nstepsRef-1)-step_size/2; %?
%             %refWt=refWt-delta*reft/refdt;
%             int_term=refdt/2*exp(0.5)+refdt*cumsum(exp(0.5*reft+refWt)) -refdt/2*exp(0.5*reft(end)+refWt(end));
%             refXt = exp(0.5*reft+refWt)./(2+int_term);                        
%             solsRef(l)=refXt(end); 
% %             figure(1);
% %             hold on;
% %             plot(reft, refWt);
% %             figure(2);
% %             hold on;
% %             plot(reft,refXt);
%           end
%         end
        
%         figure(3);
%         hist(solsRef,20);
        
        solsRef = squeeze(solRefAll(1,end,:));
        fprintf("\n mean reference \t mean SDC \t err mean \t var reference \t var SDC \t err var\n");
        fprintf("---------------------------------------------------------------------------------------------------\n");
        errMean =  abs(mean(solsRef)-mean(solsSDC));
        errVar =  abs(var(solsRef)-var(solsSDC));
        fprintf("%d \t %d \t %d \t %d \t %d \t %d\n", mean(solsRef), mean(solsSDC), errMean, var(solsRef), var(solsSDC), errVar);
        fprintf("---------------------------------------------------------------------------------------------------\n\n");
        
        %calculate weak, strong, L2-error and absolute error criterion for SDC method
        
        if strcmp(d, 'OU')
          solEx = sol(:,end);
          sumSolT2 = sumSolT2 / realIter;
          refEx= exp(-t_end);
          refEx2 = 0.5*(1+exp(-2*t_end));
          errWeak = abs(solEx-refEx);   %weak convergence, test function phi = x
          errL2 = abs(sumSolT2-refEx2); %weak convergence, test function phi = x^2
          errStrong = NaN;
          errT = NaN;
        elseif strcmp(d, 'Mattingly')          
          sumSolT2 = sumSolT2 / realIter;
          x0 = initial();
          refEx2 = x0^2*exp((2*lambda+1)*t_end) + 4*(exp((2*lambda+1)*t_end)-1)/(2*lambda+1);
          errWeak = NaN;
          errL2 = abs(sumSolT2-refEx2); %weak convergence, test function phi = x^2
          errStrong = NaN;
          errT = errL2 / refEx2;   % relative error
        else          
          errWeak = max(abs(sol-solExactRef),[],2);
          errStrong = (diffSol(:,end)./realIter);
          errL2 = (diff/realIter).^0.5;
          errT = abs(sol(:,end)-solExactRef(:,end));
        end
        


        
elseif strcmp(sde_solver,'SDC_BB_impl')
    
    parameters = [col_points, intervals, max_iter, t_begin, t_end];
    
    %computes spectral matrix (S), the micro time nodes (t) and macro time nodes (x)
    [S, t, ~] = computeSpecMat(t_begin, t_end, step_size, intervals, col_points, nodes);

        sumSol=zeros(order, nTime);
        sumSolRef=zeros(order, nTime);
        diffSol = zeros(order, nTime);
        diff = 0;
        
        compTime = 0;       
        etaCurrSteps = zeros(order, nTime);
        eta0 = zeros(order, nTime);
        xi_l = zeros(order, m-1, intervals);       
        
        for l=1:realIter
            for k=1:order
                etaMat = eta{k}(l,:);
                etaCurrSteps(k,:) = [0 etaMat(downsampleFactor:downsampleFactor:end)];
                eta0(k,:) = etaCurrSteps(k,1:end)- [0, etaCurrSteps(k,1:end-1)];
                xi_l(k,:,:) = xi{k}(l,:,downsampleFactor:downsampleFactor:end);
            end
            eta0 = sqrt(1/step_size)*eta0;
            % eta0 now has variance 1 (before: variance deltaT)
            
            if(m>1)
                %compute quadrature matrix for each expansion term
                [quadMatK_c] = quadMatKFun(t, step_size, nodes, m-1);
            else
                quadMatK_c = 0;
                %quadMatK_s = 0;
            end     
            
            %dW approximation
            deltaW = etaCurrSteps(:,2:end)- etaCurrSteps(:,1:end-1);
            
            %calls SDC_SDE_BB.m to compute SDC approximation
            tic
            [solSDC, countRhsEvaluationsThisRun] = SDC_SDE_BB_impl(parameters, initial(), S, quadMatK_c, t,step_size, nodes, beta, rhs, stochRhs, eta0, deltaW, RhsIto, order, xi_l, strInit, m, tol, J);
            elapsedTime = toc;
            compTime = compTime +elapsedTime;
            countRhsEvaluations = countRhsEvaluations + countRhsEvaluationsThisRun;
            solSDC =  solSDC(:,1:col_points-1:end);
            
            solRef = solRefAll(:,:,l);
            
            %adding values of each realization realIter
            %preliminary work for determining the expected solution
            sumSol = sumSol + solSDC;
            sumSolRef = sumSolRef + solRef;
            
            %preliminary work for computation of different error definitions
            diffSol = diffSol + abs(solSDC-solRef);
            diff = diff + (solSDC(:,end)-solRef(:,end)).^2;
            
        end %for loop realIter
        
        %calculate expected solution for each time point t
        sol= sumSol/realIter;
        solExactRef = sumSolRef/realIter;
        
        %calculate weak, strong, L2-error and absolute error criterion for SDC
        %method
        errWeak = max(abs(sol-solExactRef),[],2);
        errStrong = (diffSol(:,end)./realIter);
        errL2 = (diff/realIter).^0.5;
        errT = abs(sol(:,end)-solExactRef(:,end));
        
    
    elseif strcmp(sde_solver, 'EM')
        %% Euler Maruyama
        if (strcmp (nodes, 'equidist'))
            
            sumSol=zeros(order, nTime);
            sumSolRef=zeros(order, nTime);
            diffSol = zeros(order, nTime);
            diff = 0;
            compTime = 0;
            etaCurrSteps = zeros(order, nTime);
            
            %loop over #realizations
            for l=1:realIter
                %choose Brownian motion & compute its increments
                for k=1:order
                    etaMat = eta{k}(l,:);
                    etaCurrSteps(k,:) = [0 etaMat(downsampleFactor:downsampleFactor:end)];
                end
                deltaW = etaCurrSteps(:,2:end) - etaCurrSteps(:,1:end-1);
                
                %calls EulerMaruyama.m
                tic
                [ sol, countRhsEvaluationsThisRun, ~ ] = EulerMaruyama( initial, step_size, t_begin, t_end, RhsIto, stochRhs, deltaW, order);
                elapsedTime = toc;
                %compute time for EM approximation
                compTime = compTime +elapsedTime;
                countRhsEvaluations = countRhsEvaluations + countRhsEvaluationsThisRun;
                
                sumSol = sumSol + sol;
                %exact solution at macro time steps
                solRef = solRefAll(:,:,l);
                sumSolRef = sumSolRef + solRef;
                
                %basis to compute strong convergence error
                diffSol = diffSol + abs(sol-solRef);
                diff = diff + (sol(:,end)-solRef(:,end)).^2;
                
            end %realIter
            
            %expectation of the approximation and solution
            solExactRef = sumSolRef/realIter;
            sol = sumSol/realIter;
            
            %calculate weak, strong, L2-error and absolute error criterion for
            %approximation method
            errWeak = max(abs(sol-solExactRef),[],2);
            errStrong = diffSol(:,end)./realIter;
            errL2 = (diff/realIter).^0.5;
            errT = abs(sol(:,end)-solExactRef(:,end));
        else
            fprintf('Euler Maruyama: choose equidistant time grid');
            return;
        end
        
    elseif strcmp(sde_solver, 'Mil')
        %% Milstein method
        if (strcmp (nodes, 'equidist'))
            
            sumSol=zeros(order, nTime);
            sumSolRef=zeros(order, nTime);
            diffSol = zeros(order, nTime);
            diff = 0;
            compTime = 0;
            etaCurrSteps = zeros(order, nTime);
            
            for l=1:realIter
                %choose Brownian motion and compute its increments
                for k=1:order
                    etaMat = eta{k}(l,:);
                    etaCurrSteps(k,:) = [0 etaMat(downsampleFactor:downsampleFactor:end)];
                end
                deltaW = etaCurrSteps(:,2:end) - etaCurrSteps(:,1:end-1);
                
                %calls Milstein.m
                %notice that in multidim case, Milstein uses the Stratonovich
                %form not the Ito one
                if (order >1)
                    tic
                    [ sol, countRhsEvaluationsThisRun, ~ ] = Milstein( initial, step_size, t_begin, t_end, beta, rhs, stochRhs, J, deltaW, order);
                    elapsedTime = toc;
                else
                    tic
                    [ sol, countRhsEvaluationsThisRun, ~ ] = Milstein( initial, step_size, t_begin, t_end, beta, RhsIto, stochRhs, J,deltaW, order);
                    elapsedTime = toc;
                end
                compTime = compTime +elapsedTime;
                countRhsEvaluations = countRhsEvaluations + countRhsEvaluationsThisRun;
                
                sumSol = sumSol + sol;
                %exact solution at macro time steps
                solRef = solRefAll(:,:,l);
                sumSolRef = sumSolRef + solRef;
                
                %for strong convergence error
                diffSol = diffSol + abs(sol-solRef);
                diff = diff + (sol(:,end)-solRef(:,end)).^2;
                
            end %realIter
            
            %expectation of solution and approximation
            solExactRef = sumSolRef/realIter;
            sol = sumSol/realIter;
            
            %calculate weak, strong, L2-error and absolute error criterion for
            %approximation method
            errWeak = max(abs(sol-solExactRef),[],2);
            errStrong = diffSol(:,end)./realIter;
            errL2 = (diff/realIter).^0.5;
            errT = abs(sol(:,end)-solExactRef(:,end));
        else
            fprintf('Milstein: choose equidistant time grid');
            return;
        end
    elseif strcmp(sde_solver, 'implMilDampNM') %%%WORKS ONLY FOR TP2!!!
        %% implicit Milstein method
        if (strcmp (nodes, 'equidist'))
            
            sumSol=zeros(order, nTime);
            sumSolRef=zeros(order, nTime);
            diffSol = zeros(order, nTime);
            diff = 0;
            compTime = 0;
            etaCurrSteps = zeros(order, nTime);
            
            for l=1:realIter
                %choose Brownian motion and compute its increments
                for k=1:order
                    etaMat = eta{k}(l,:);
                    etaCurrSteps(k,:) = [0 etaMat(downsampleFactor:downsampleFactor:end)];
                end
                deltaW = etaCurrSteps(:,2:end) - etaCurrSteps(:,1:end-1);
                
                %calls Milstein.m
                %notice that in multidim case, Milstein uses the Stratonovich
                %form not the Ito one
                if (order >1)
                    tic
                    [ sol, countRhsEvaluationsThisRun, ~ ] = ImplMilsteinDampNM( initial, step_size, t_begin, t_end, beta, rhs, stochRhs, J, deltaW, order);
                    elapsedTime = toc;
                else
                    tic
                    [ sol, countRhsEvaluationsThisRun, ~ ] = ImplMilsteinDampNM( initial, step_size, t_begin, t_end, beta, RhsIto, stochRhs, J,deltaW, order);
                    elapsedTime = toc;
                end
                compTime = compTime +elapsedTime;
                countRhsEvaluations = countRhsEvaluations + countRhsEvaluationsThisRun;
                
                sumSol = sumSol + sol;
                %exact solution at macro time steps
                solRef = solRefAll(:,:,l);
                sumSolRef = sumSolRef + solRef;
                
                %for strong convergence error
                diffSol = diffSol + abs(sol-solRef);
                diff = diff + (sol(:,end)-solRef(:,end)).^2;
                
            end %realIter
            
            %expectation of solution and approximation
            solExactRef = sumSolRef/realIter;
            sol = sumSol/realIter;
            
            %calculate weak, strong, L2-error and absolute error criterion for
            %approximation method
            errWeak = max(abs(sol-solExactRef),[],2);
            errStrong = diffSol(:,end)./realIter;
            errL2 = (diff/realIter).^0.5;
            errT = abs(sol(:,end)-solExactRef(:,end));
        else
            fprintf('impl Milstein: choose equidistant time grid');
            return;
        end
        
    else
        %% implicit Euler Maruyama with damped Newton's method
        % 'implEMDampNM'
        if (strcmp (nodes, 'equidist'))
            
            sumSol=zeros(order, nTime);
            sumSolRef=zeros(order, nTime);
            diffSol = zeros(order, nTime);
            diff = 0;
            compTime = 0;
            etaCurrSteps = zeros(order, nTime);
            
            %loop over #realizations
            for l=1:realIter
                %choose Brownian motion & compute its increments
                for k=1:order
                    etaMat = eta{k}(l,:);
                    etaCurrSteps(k,:) = [0 etaMat(downsampleFactor:downsampleFactor:end)];
                end
                deltaW = etaCurrSteps(:,2:end) - etaCurrSteps(:,1:end-1);
                
                %calls EulerMaruyama.m
                tic
                [ sol, countRhsEvaluationsThisRun, ~ ] = ImplEulerMaruyamaDampNM( initial, step_size, t_begin, t_end, RhsIto, stochRhs, deltaW, order);
                elapsedTime = toc;
                %compute time for EM approximation
                compTime = compTime +elapsedTime;
                countRhsEvaluations = countRhsEvaluations + countRhsEvaluationsThisRun;
                
                sumSol = sumSol + sol;
                %exact solution at macro time steps
                solRef = solRefAll(:,:,l);
                sumSolRef = sumSolRef + solRef;
                
                %basis to compute strong convergence error
                diffSol = diffSol + abs(sol-solRef);
                diff = diff + (sol(:,end)-solRef(:,end)).^2;
                
            end %realIter
            
            %expectation of the approximation and solution
            solExactRef = sumSolRef/realIter;
            sol = sumSol/realIter;
            
            %calculate weak, strong, L2-error and absolute error criterion for
            %approximation method
            errWeak = max(abs(sol-solExactRef),[],2);
            errStrong = diffSol(:,end)./realIter;
            errL2 = (diff/realIter).^0.5;
            errT = abs(sol(:,end)-solExactRef(:,end));
        else
            fprintf('Euler Maruyama: choose equidistant time grid');
            return;
        end
        
end %choice of method
    
    if (plot_Error == true)
        plot_err(step_size, errWeak, errStrong, errL2, errT);
    end
    if (plot_Sol == true)
        plot_sol(time, t_begin, t_end, intervals, sol, solExactRef);
    end
    
    %fprintf('\n end of sdeMethod.m \n ');
    
end %function
