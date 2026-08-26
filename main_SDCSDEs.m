%% preperation for applying the chosen method
%generates data files
%by Lisa Fischer
function [] = main_SDCSDEs(d, sde_solver, nodes, strInit, colpoints, maxIter, steps, NNfinest, realIter, mBB, plot_Sol, plot_Error, SBB)

fprintf('START PROGRAM using the solver: %s \n', sde_solver)

seed = 2348;
rng(seed)

%%setting
%time interval & step_size
t_begin = 0;
t_end = 1;

%stopping criterion for SDC method (additionally to #sweeps)
tol = 1e-10;

%Definition of the SDE
if strcmp(d, 'exp')
    lambda = 0.5;
    beta = 1.0;
    [rhs, stochRhs, J, RhsIto, exact] = problem1(lambda, beta);
    initial = 1.0;
    order = 1;
elseif strcmp(d, 'TP2')
    lambda = 1;
    beta = 1;
    [rhs, stochRhs, J, RhsIto, exact] = problem2(lambda, beta);
    initial = 0.0;
    order = 1;
elseif strcmp(d, 'TP3')
    lambda = 1;
    beta = 1;
    [rhs, stochRhs, J, RhsIto, exact] = problem3(lambda, beta);
    initial = 0.5;
    order = 1;
    t_end = steps(1); % for doing one timestep only
elseif strcmp(d, 'TP4')
    lambda = [1.5 -0.85; 1.275 -0.625];
    B1 = [0.9 -0.2; 0.3 0.4];
    B2 = [2.1 -1.2; 1.8 -0.9];
    beta = {B1, B2};
    [rhs, stochRhs, J, RhsIto, exact] = problem4(lambda, B1, B2);
    initial = [1; 0];
    order = 2;
elseif strcmp(d, 'OU')
    lambda = -0.5;
    beta = 1.0;
    [rhs, stochRhs, J, RhsIto, exact] = problemOU(lambda, beta);
    initial = 0.5;
    order = 1;
elseif strcmp(d, 'Mattingly')
    lambda = -2.0;
    beta = 1.0;
    [rhs, stochRhs, J, RhsIto, exact] = problemMattingly(lambda, beta);
    initial = 2.0;
    order = 1;
    t_end = 2;
else
    lambda =0;
    beta = 1;
    [rhs, stochRhs, J, RhsIto, exact] = problem(lambda, beta);
    initial = 0.0;
    order = 1;
end

%time grid
nSteps = length(steps);
fprintf('finest time grid: %d', NNfinest);
nRndVar = NNfinest;
eta = cell(1,order);
xi = cell(1,order);
%Brownian motion W(t) valid on [0,1] without W(0)
for k=1:order
    %eta{k} = cumsum(randn(realIter,nRndVar),2).*1./sqrt(nRndVar); % TAKE CARE: only valid for t_end=1!
    eta{k} = cumsum(randn(realIter,nRndVar),2).*sqrt(t_end)./sqrt(nRndVar); 
end

%initialization for error output files
errW = zeros(order,nSteps);
errStr = zeros(order,nSteps);
errorL2 = zeros(order,nSteps);
errorT = zeros(order,nSteps);

%global solRefFinestAll

if strcmp(d, 'OU') || strcmp(d, 'Mattingly')
  solRefFinestAll = 0;
else
  solRefFinestAll = zeros(order, NNfinest+1, realIter);
end 
etaFin = zeros(order, NNfinest+1);
etaFinest = zeros(order, NNfinest+1);

for p=1:length(mBB)
    %loops over different many Karhunen-Loeve expansion terms
    m = mBB(p);
    if (m >1)
        %draw samples for more than one Karhunen-Loeve expansion term
        for k=1:order
            xi{k} = randn(realIter,m-1, NNfinest);
        end
    else
        for k=1:order
            xi{k} = zeros(realIter,0, NNfinest);
        end
    end
    xiFinest = zeros(order, m-1, NNfinest);
    
    
    %save output in folder with name of #Karhunen-Loeve expansion terms,
    %#realizations and finest time grid
    subfolder = sprintf('data_bb_%d_%d_%d', m, realIter, NNfinest);
    if SBB == true
      subfolder = sprintf('data_bb_SBB_%d_%d_%d', m, realIter, NNfinest);
    end 
    mkdir(subfolder);
    
    if strcmp(d, 'OU') || strcmp(d, 'Mattingly')
      solRefFinestAll = 0;
      fprintf(1, '\n problem: %s \t no ref sol required' , d);
%    elseif strcmp(d, 'TP3')      
%      allSolRef = cell(1,length(steps));
%      for n=1:length(steps)
%        allSolRef{n} = zeros(order, round((t_end-t_begin)/steps(n))+1, realIter);
%      end
%      fprintf(1, '\n problem: %s \t ref sol will be computed on current sample path, not on finest' , d);
%      
% %       % check bridge ODE with Matlab solver
%       for l=1:realIter
%         stepsFinest = t_end/NNfinest;
%         for k=1:order
%           etaMat = eta{k}(1,:);%eta{k}(l,:); %always use one to have all processes go through the same nodes
%           etaFin(k,:) = [0 etaMat];
%           etaFinest(k,:) = etaFin(k,1:end)- [0, etaFin(k,1:end-1)];
%           xiFinest(k,:,:) = xi{k}(l,:,:);
%         end
%         etaFinest = sqrt(1/stepsFinest)*etaFinest; % t_end/stepsFinest?
%         if l==1
%           fprintf(1, '\n etaFinest= %d' , etaFinest(1,end));
%         end
%         k=1;
%         bM(:,k) = brownianBridge(etaFinest(:,k+1), stepsFinest , stepsFinest, xiFinest);
%         for k=2:NNfinest
%            bM(:,k) =  brownianBridge(etaFinest(:,k+1), stepsFinest , stepsFinest, xiFinest)+bM(:,k-1);
%         end
%         bbM(:,:) = [zeros(order,1) bM(:,:)];
%         y0 = initial;
%         tsolAppro = t_begin;
%         solAppro = initial;
%         for k=1:NNfinest
%           [tsol, solA] = ode45(@(s,x)SBB_ODE_TP3(s,x, etaFinest(:,k+1), stepsFinest, xiFinest(:,:,k), lambda, beta),[t_begin+(k-1)*stepsFinest t_begin+k*stepsFinest], y0');
%           y0 = solA(end,:);
%           tsolAppro = [tsolAppro tsol(end)'];
%            solAppro = [solAppro solA(end)'];
%         end
%         solRefFinestAll(:,:,l) = solAppro;
% 
% %         if mod(l,100000) == 0
% %           figure(10);
% %           hold on ;
% %           plot(1:size(etaFinest,2), squeeze(etaFinest(1,:)));
% %         end
%       end
%       sv = squeeze(solRefFinestAll(1,end,:));
%       mean(sv)
%       var(sv)
%       return;
    else
    
    %sol of the smooth Brownian Bridge ODE
    % on finest time grid
    for l=1:realIter
        stepsFinest = t_end/NNfinest;
        for k=1:order
            etaMat = eta{k}(l,:);
            etaFin(k,:) = [0 etaMat];
            etaFinest(k,:) = etaFin(k,1:end)- [0, etaFin(k,1:end-1)];
            xiFinest(k,:,:) = xi{k}(l,:,:);
        end
        etaFinest = sqrt(1/stepsFinest)*etaFinest; % t_end/stepsFinest?
        
        if SBB == true
        k=1;
        bM(:,k) = brownianBridge(etaFinest(:,k+1), stepsFinest , stepsFinest, xiFinest);
        for k=2:NNfinest
            bM(:,k) =  brownianBridge(etaFinest(:,k+1), stepsFinest , stepsFinest, xiFinest)+bM(:,k-1);
        end
        
        bbM(:,:) = [zeros(order,1) bM(:,:)];
        solRefFinestAll(:,:,l) = exact(lambda, beta, t_begin, t_end, stepsFinest, initial, bbM);
              %ref sol using matlab ode solver  
%                 y0 = initial;
%                 tsolAppro = t_begin;
%                 solAppro = initial;
%                 for k=1:NNfinest
%                     %[tsol, solA] = ode45(@(s,x)SBB_ODE_TP3(s,x, etaFinest(:,k+1), stepsFinest, xiFinest(:,:,k), lambda, beta),[t_begin+(k-1)*stepsFinest t_begin+k*stepsFinest], y0');
%                     [tsol, solA] = ode45(@(s,x)SBB_ODE_exp(s,x, etaFinest(:,k+1), stepsFinest, xiFinest(:,:,k), lambda, beta, k),[t_begin+(k-1)*stepsFinest t_begin+k*stepsFinest], y0');
%                     y0 = solA(end,:);
%                     tsolAppro = [tsolAppro tsol(end)'];
%                     solAppro = [solAppro solA(end)'];
%                 end
%                 solRefFinestAll(:,:,l) = solAppro;
        else
            solRefFinestAll(:,:,l) = exact(lambda, beta, t_begin, t_end, stepsFinest, initial, etaFin);
        end
        
    end
    end
    
    
    %loop over all collocation points
    for k=1:length(colpoints)
        
        data = cell(1,order);
        for s=1:order
            data{s}=zeros(nSteps,9);
        end
        
        col_points = colpoints(k);
        max_iter = maxIter(k);
        fprintf(1, '\n problem: %s \t nodes: %s \t time interval: [%d, %d] \t col points: %d \t max iter: %d \t m= %d' , d, nodes, t_begin, t_end, col_points, max_iter, m);
        
        %loop over all time steps
        for n=1:length(steps)                        
            step_size = steps(n);
            intervals = round((t_end-t_begin)/step_size);
            fprintf(1, '\nstep size: %d\t intervals: %d', steps(n), intervals);
            
            %global solRef
            if strcmp(d, 'OU') || strcmp(d, 'Mattingly')
              solRef = 0;
            elseif strcmp(d, 'TP3')
              
%               if k == 1            
%                 etaCurrSteps = zeros(order,intervals+1);
%                 for l=1:realIter
%                   for o=1:order
%                     etaMat = eta{o}(l,:);
%                     etaCurrSteps(o,:) = [0 etaMat(NNfinest/intervals:NNfinest/intervals:end)];                   
%                   end         
%                   allSolRef{n}( :, :, l) = exact(lambda, beta, t_begin, t_end, step_size, initial, etaCurrSteps);
%                 end
%               end
              solRef = zeros(order, intervals+1, realIter);
%               solRef(:,:,:) = allSolRef{n}(:,:,:);
            else
              solRef = zeros(order, intervals+1, realIter);
              solRef(:,:,:) = solRefFinestAll(:,1:NNfinest/intervals:end,:);                
            end
            
            nameDat = sprintf('%s_%s_%d_%s_%s_%d%d_%d_%d_%d_%d', sde_solver, d, col_points, nodes(1:2), strInit(1:6), t_begin, t_end, realIter, NNfinest, m, max_iter);
            
            %function call to compute EM, Milstein or SDC approximation
            [errWeak, errStrong, errL2, errT, compTime, countRhsEvaluations] = sdeMethod(sde_solver, tol, NNfinest, step_size, t_begin, t_end, intervals, col_points, max_iter, realIter, initial, nodes, lambda, beta, rhs, stochRhs, RhsIto, J, order, eta, exact, d, xi, plot_Error, plot_Sol, strInit,m, solRef);
            
            %preparation for output saving
            errW(:,n) = errWeak;
            errStr(:,n) = errStrong;
            errorL2(:,n) = errL2;
            errorT(:,n) = errT;
            
            for l=1:order
                data{1,l}(n, :) = [step_size, col_points, max_iter, errWeak(l), errStrong(l), errL2(l), errT(l), countRhsEvaluations, compTime];
            end
            
            %end loop over time steps if close to machine precision
            if (errStrong(1) < 10^-12 && errWeak(1) < 10^-12)
                n = length(steps)+1;
                break;
            end
        end %for loop step size
        
        %errW, errStr, errorL2, errorT, errWapp, errStrApp, errorL2app, errorTapp,
        %compute slope as an indicator for convergence order
        %does not work if loop over time steps ends before length(step_size)
%         slopeStrong = (log(errStr(1))-log(errStr(end)))/(log(steps(1))-log(steps(end)));
%         slopeWeak = (log(errW(1))-log(errW(end)))/(log(steps(1))-log(steps(end)));
%         slopeL2 = (log(errorL2(1))-log(errorL2(end)))/(log(steps(1))-log(steps(end)));
%         slopeT = (log(errorT(1))-log(errorT(end)))/(log(steps(1))-log(steps(end)));
        
        eocStr = zeros(size(errStr));
        eocWeak = zeros(size(errW));
        eocT = zeros(size(errorT));
        eocL2 = zeros(size(errorL2));
        for i=1:size(errStr,2)-1
          eocStr(:,i+1) = (log(errStr(:,i+1))-log(errStr(:,i))) / (log(steps(i+1))-log(steps(i)));
          eocWeak(:,i+1) = (log(errW(:,i+1))-log(errW(:,i))) / (log(steps(i+1))-log(steps(i)));
          eocT(:,i+1) = (log(errorT(:,i+1))-log(errorT(:,i))) / (log(steps(i+1))-log(steps(i)));
          eocL2(:,i+1) = (log(errorL2(:,i+1))-log(errorL2(:,i))) / (log(steps(i+1))-log(steps(i)));
        end
        
        fprintf( '\n convergence order: \t col points: %d \t max iter: %d \n', col_points, max_iter);
        %fprintf( 'slope: errStrong %f \t, errWeak %f \t, errL2 %f \t, errT %f \n', slopeStrong, slopeWeak, slopeL2, slopeT)
        fprintf( '\n step size \t strong err \t eoc \t weak err \t eoc \t L2 err \t eoc \t T err \t eoc\n');
        fprintf( '----------------------------------------------------------------------------------------------------------------------------------\n');
        for i=1:length(steps)
          fprintf( '%1.6e \t %1.6e \t %1.2f \t %1.6e \t %1.2f \t %1.6e \t %1.2f \t %1.6e \t %1.2f  \n', steps(i), errStr(1,i), eocStr(i), errW(1,i), eocWeak(i), errorL2(1,i), eocL2(i), errorT(1,i), eocT(i));
        end
        
        fprintf( '\n compTime %f\n' , compTime);
        
        %save error data in file
        for q=1:order
            if (order >1)
                filename = sprintf('Err_compEff_X%d_%s.dat', q,nameDat);
            else
                filename = sprintf('Err_compEff_%s.dat', nameDat);
            end
            fullFileName = fullfile(subfolder, filename);
            fileID = fopen(fullFileName,'w');
            fprintf(fileID, 'stepSize \t colpoints \t maxIter \t errW \t errStr \t errorL2 \t errT \t fctEval \t compTime \n');
            fclose(fileID);
            
            %dat = cell2mat(data{s})
            fileID = fopen(fullFileName,'a+');
            fprintf(fileID, '%f \t %d \t %d \t %d \t %e \t %e \t %e \t %d \t %e \n', data{q}');
            fclose(fileID);
        end
        
    end %for loop colpoints
    
end %for loop mBB

end