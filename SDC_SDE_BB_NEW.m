%% sequential Spectral Deferred Correction Method
% using the explicit Euler scheme

% parameters - number of collocation points, intervals, max sweeps
% initialValue - initial value for SDE
% S - quadrature matrix
% quadMatK_c - quadrature matrix for the cosines times monomial
% t - time intervals (Wiener process sampled at these time steps)
% deltaT - time step size
% nodes - string specifying type of collocation nodes
% beta - coefficients for stochastic part (length = # independent Wiener processes)
% rhs_eval - count evaluations of rhs (drift and diffusion)
% stochRhs - function handle for the diffusion part
% eta - random variables specifying the bridge
% deltaW - increments of the Wiener process
% RhsIto - drift term in Ito formulation (for some initialization method)
% order - dimension of the system (number of unknowns)
% xi - random variables specifying the bridge
% strInit - string specifying type of intialization
% nBridgeTerms - number of terms in the bridge
% tol - stopping tolerance: stop sweeps if correction norm < tol

function [ sol, countRhsEvaluations ] = SDC_SDE_BB_NEW( parameters, initialValue, S, quadMatK_c, t ,deltaT, nodes , beta, rhs_eval, stochRhs, eta, deltaW, RhsIto, order, xi, strInit, nBridgeTerms, tol)

col_points = parameters(1);
intervals = parameters(2);
max_sweeps = parameters(3);


%initials & placeholders needed for implementation
if strcmp (nodes, 'radau')
    points = col_points+1;
elseif strcmp (nodes, 'legendre')
    points = col_points+2;
else
    points = col_points;
end

% initial values needed for SDC method
d = zeros(order, points);
sol = zeros(order, (points-1)*intervals+1);

% initial value of SDE
y0 = initialValue;

correctionNorm=-1;
countRhsEvaluations = 0;

t_begin = t(1);
%for each sub time interval
for i=1:intervals
    %fprintf('\n %d interval of %d \n', i, intervals)
    
    %draw samples for more than one Karhunen-Loeve expansion term
    %xi = randn(order,m-1);
    
    int_begin = 1 + (i-1)* (points-1);
    int_end = int_begin + points-1   ;
    
    %compute derivative B0 Brownian Bridge term -> per definition
    %eta/sqrt(deltaT), with eta ~N(0,1)
    db0 = eta(:,i+1)/sqrt(deltaT); 
    
    t_currInt = col_nodes(t_begin,t_begin+deltaT,col_points,nodes);
    
    %choose the different possible initializations
    if strcmp(strInit, 'constInit')
        phi = zeros(order, points);
        for n=1:order
            phi(n, :) = ones(1,points)*y0(n);
        end
        
        if(nBridgeTerms>1)
            %dbBridge = dbrownianBridge(eta(:,i+1), deltaT, 0, xi(:,:,i));
            dbBridge = dbrownianBridge(deltaT, t_currInt(1), xi(:,:,i));
        else
            dbBridge = zeros(order,1);
        end
        %Brownian bridge
        countRhsEvaluations = countRhsEvaluations + 1;
        
    elseif strcmp(strInit, 'eulerSBBInit')
        % set initial value for current time interval
        %phi(:,1) = y0;
        rhs_phi = rhs_eval(y0);
        phi(:,1) = y0;
        
        deltaTau = t(2:end)-t(1:end-1);
        %rhs_phi fct evals
        countRhsEvaluations = countRhsEvaluations + 1;
        
        for k=2:points
            if(nBridgeTerms>1)
                dbBridge = dbrownianBridge(deltaT, t_currInt(k-1), xi(:,:,i));
            else
                dbBridge = zeros(order,1);
            end
            %Brownian bridge
            countRhsEvaluations = countRhsEvaluations + 1;
            
            %stoch part
            stoch_part = zeros(order,1);
            for n=1:size(beta,2)
                stoch_part = stoch_part + stochRhs{n}(phi(:,k-1)).*(db0(n) + dbBridge(n));
                %rhs stoch fct evals
                countRhsEvaluations = countRhsEvaluations + 1;
            end
            %using one Euler-SBB step
            phi(:,k) = phi(:,k-1) + deltaTau(k-1)*rhs_phi + deltaTau(k-1)*stoch_part;
            rhs_phi = rhs_eval(phi(:,k));
            %rhs_phi fct evals
            countRhsEvaluations = countRhsEvaluations + 1;
        end
    else %EM Ito SDE appl linear interpolation
        % set initial value for current time interval
        
        phi= y0;
        rhs_phi = RhsIto(phi);
        %rhs_phi fct evals
        countRhsEvaluations = countRhsEvaluations + 1;
        
        %stoch part
        stoch_part = zeros(order,1);
        for n=1:size(beta,2)
            stoch_part = stoch_part + stochRhs{n}(phi(:,1))*deltaW(n,i);
            %rhs_stoch fct evals
            countRhsEvaluations = countRhsEvaluations + 1;
        end
        
        for k=2:points
            %using one Euler-Maruyama step and linear interpolation for all
            %collocation points
            phi(:,k) = phi(:,1) + (rhs_phi*deltaT + stoch_part)*t(k)/deltaT;
            
        end
        
        if(nBridgeTerms>1)
            %dbBridge = dbrownianBridge(eta(:,i+1), deltaT, 0, xi(:,:,i));
            dbBridge = dbrownianBridge(deltaT, t_currInt(1), xi(:,:,i));
        else
            dbBridge = zeros(order,1);
        end
        %Brownian bridge
        countRhsEvaluations = countRhsEvaluations + 1;
    end %choice of initialization
    
    
    %maximum number of sweeps
    for j=1:max_sweeps
        
        % SDC scheme for integration over all collocation points
        for p=2:points
            
            % deterministic part - 2 fct evals
            rhs = rhs_eval(phi);
            rhs_integrate = (S(p-1,:)*rhs')';
            rhs_err = rhs_eval(phi(:,p-1)+d(:, p-1));
            rhs_phi_prev = rhs(:,p-1);
            countRhsEvaluations = countRhsEvaluations + 2;
            
            %stoch part
            stoch_rhs_integrate_b0 = zeros(order, 1);
            stoch_rhs_integrate_bm = zeros(order, 1);
            stoch_rhs_diff = zeros(order, 1);
            for n=1:size(beta,2)
                %terms for stoch spectral integration
                stoch_rhs = stochRhs{n}(  phi());
                stoch_rhs_err = stochRhs{n}(  phi(:,p-1)+d(:, p-1) );
                stoch_rhs_prev = stoch_rhs(:,p-1);
                countRhsEvaluations = countRhsEvaluations + 2;
                
                stoch_rhs_integrate_b0 = stoch_rhs_integrate_b0 + db0(n)*(S(p-1,:)*stoch_rhs')';
                if (nBridgeTerms>1)
                    % loop over expansion terms: sum_1^m {xi_k * Q_k}
                    % include random variables xi as in the Karhunen-Loeve expansion
                    Q = zeros(order,size(quadMatK_c, 2));
                    for l=1:nBridgeTerms-1
                        Q = Q + quadMatK_c(p-1,:, l).*xi(:,l,i);
                    end
                    
                    Q= sqrt(2/deltaT)*Q;
                    for y=1:order
                    stoch_rhs_integrate_bm(y) = (Q(y,:)*stoch_rhs(y,:)')';
                    end
                    dbBridge = dbrownianBridge(deltaT, t_currInt(p-1), xi(:,:,i));
                else
                    stoch_rhs_integrate_bm = zeros(order,1);
                end
                countRhsEvaluations = countRhsEvaluations + 1;
                
                stoch_rhs_diff = stoch_rhs_diff + (db0(n) + dbBridge(n))*(stoch_rhs_err - stoch_rhs_prev);
            end %end loop system's order
            
            %correction scheme
            d(:,p) = d(:, p-1) + (t(:,p) - t(:,p-1)) * ( rhs_err - rhs_phi_prev) + ...
                (t(:,p) - t(:,p-1)) * stoch_rhs_diff...
                + rhs_integrate - phi(:,p) + phi(:, p-1) ...
                + stoch_rhs_integrate_b0 + stoch_rhs_integrate_bm;
            
        end
        
        %solution update
        phi = phi + d;
        
        correctionNorm = norm(d);
        
        if (correctionNorm < tol)
            break;
        end
        
    end %iterations
    if correctionNorm > tol
        correctionNorm
    end
    y0 = phi(:,end);
    sol(:,int_begin:int_end) = phi;
    d = zeros(order, points);
    t_begin = t_currInt(end);
    
end %intervals


end %function
