function [] = main_realIter_moments(d, sde_solver, nodes, strInit, colpoints, maxIter, steps, NNfinest, realIter, mBB, tolMoments)

seed = 2348;
rng(seed)

%time interval & step_size
t_begin = 0;
t_end = 1;

tol = 1e-12;

%Definition of the SDE
if strcmp(d, 'exp')
    lambda = 0.5;
    beta = 1.0;
    [rhs, stochRhs, J, RhsIto, exact] = problem1(lambda, beta);
    initial = 1.0;
    order = 1;
elseif strcmp(d, 'TP2')
    lambda = 0;
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
else
    lambda = [1.5 -0.85; 1.275 -0.625];
    B1 = [0.9 -0.2; 0.3 0.4];
    B2 = [2.1 -1.2; 1.8 -0.9];
    beta = {B1, B2};
    [rhs, stochRhs, J, RhsIto, exact] = problem4(lambda, B1, B2);
    initial = [1; 0];
    order = 2;
end

%time grid
nSteps = length(steps);
fprintf('finest time grid: %d', NNfinest);
nRndVar = NNfinest;
eta = cell(1,order);
xi = cell(1,order);
%Brownian motion W(t) valid on [0,1] without W(0)
for k=1:order
    eta{k} = cumsum(randn(realIter,nRndVar),2).*1./sqrt(nRndVar); % TAKE CARE: only valid for t_end=1!
end


for q=1:length(tolMoments)
    tolMom = tolMoments(q);
for p=1:length(mBB)
    %loops over different many Karhunen-Loeve expansion terms
    m = mBB(p);
    if (m >1)
        %draw samples for more than one Karhunen-Loeve expansion term
        for k=1:order
            xi{k} = randn(realIter,m-1);
        end
    else
        for k=1:order
           xi{k} = zeros(realIter,1);
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
              
            %function call to compute realIter for EM, Milstein or SDC
            %approximation reaching tol_moments
            realIterMoments(sde_solver, tol, NNfinest, step_size, t_begin, t_end, intervals, col_points, max_iter, realIter, initial, nodes, lambda, beta, rhs, stochRhs, RhsIto, J, order, eta, exact, d, xi, strInit,m, tolMom);
           
        end %for loop step size
              
    end %for loop colpoints
    
end %for loop mBB
end %tolMoments

end