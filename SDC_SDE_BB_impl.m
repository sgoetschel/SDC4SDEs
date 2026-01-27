%% sequential Spectral Deferred Correction Method
% using the implicit Euler-Maruyama scheme (difussive term) and explicit scheme for drift
% term

function [ sol, countRhsEvaluations ] = SDC_SDE_BB_impl( parameters, initial, S, quadMatK_c, t ,deltaT, nodes , beta, rhs_eval, stochRhs, eta, deltaW, RhsIto, order, xi_l, strInit, m, tol, J)

col_points = parameters(1);
intervals = parameters(2);
max_iter = parameters(3);


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
y0 = initial;

correctionNorm=-1;
countRhsEvaluations = 0;

t_begin = t(1);
%for each sub time interval
for i=1:intervals
    
    %draw samples for more than one Karhunen-Loeve expansion term
    %xi = randn(order,m-1);
    xi = xi_l(:,:,i);
    
    int_begin = 1 + (i-1)* (points-1);
    int_end = int_begin + points-1   ;
    
    %compute derivative B0 Brownian Bridge term
    db0 = eta(:,i+1)./sqrt(deltaT);
    
    t_currInt = col_nodes(t_begin,t_begin+deltaT,col_points,nodes);
    
    %choose the different possible initializations
    if strcmp(strInit, 'constInit')
        phi = zeros(order, points);
        for n=1:order
            phi(n, :) = ones(1,points)*y0(n);
        end
        
        if(m>1)
            dbBridge = dbrownianBridge(deltaT, 0, xi);
        else
            dbBridge = zeros(order,1);
        end
        %Brownian bridge
        countRhsEvaluations = countRhsEvaluations + 1;
        
    elseif strcmp(strInit, 'eulerSBBInit')
        % set initial value for current time interval
        rhs_phi = rhs_eval(y0);
        phi(:,1) = y0;
        
        %rhs_phi fct evals
        countRhsEvaluations = countRhsEvaluations + 1;
        
        deltaTau = t(2:end)-t(1:end-1);
        
        %Brownian bridge
        for k=2:points
            if(m>1)
                dbBridge = dbrownianBridge(deltaT, t_currInt(k-1), xi);
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
        %%%WORKS ONLY FOR TP2!!!
    elseif strcmp(strInit, 'implDampNM_eMItoSDE')
        %implicit EM Ito SDE appl linear interpolation and damped Newton's method
        % set initial value for current time interval
        rhs_phi = RhsIto(y0);
        phi(:,1) = y0;
        %rhs_phi fct evals
        countRhsEvaluations = countRhsEvaluations + 1;
        
        %stoch part
        stoch_part = zeros(order,1);
        for n=1:size(beta,2)
            stoch_part = stoch_part + stochRhs{n}(phi(:,1))*deltaW(n,i);
            %rhs_stoch fct evals
            countRhsEvaluations = countRhsEvaluations + 1;
        end
        
        %damped Newton's method for implicit Euler-Maruyama method
        damping = 1;
        phi_end = y0;
        delta_phi = -(phi_end - rhs_phi*deltaT-y0-stoch_part)/(1-(3*phi_end*phi_end-1)*deltaT);
        niter = 0;
        while abs(delta_phi)>1e-5 && niter < 10
            phi_end_new = phi_end + damping*delta_phi;
            rhs_phi = RhsIto(phi_end_new);
            countRhsEvaluations = countRhsEvaluations + 1;
            delta_phi_s = -(phi_end_new - rhs_phi*deltaT-y0-stoch_part)/(1-(3*phi_end*phi_end-1)*deltaT);
            while abs(delta_phi_s)^2 > (1-damping/2)*abs(delta_phi)^2 && damping > 1e-6
                damping = 0.5*damping;
                phi_end_new = phi_end + damping*delta_phi;
                rhs_phi = RhsIto(phi_end_new);
                countRhsEvaluations = countRhsEvaluations + 1;
                delta_phi_s = -(phi_end_new - rhs_phi*deltaT-y0-stoch_part)/(1-(3*phi_end*phi_end-1)*deltaT);
            end
            damping = min(1,2*damping);
            phi_end = phi_end_new;
            delta_phi = -(phi_end - rhs_phi*deltaT-y0-stoch_part)/(1-(3*phi_end*phi_end-1)*deltaT);
            niter = niter+1;
        end
        
        %linear interpolation for all collocation points
        for k=2:points
            phi(:,k) = y0 + (phi_end-y0) *t(k)/deltaT;
        end
        
        if(m>1)
            dbBridge = dbrownianBridge(deltaT, 0, xi);
        else
            dbBridge = zeros(order,1);
        end
        %Brownian bridge
        countRhsEvaluations = countRhsEvaluations + 1;
        
    else %EM Ito SDE appl linear interpolation
        % set initial value for current time interval
        rhs_phi = RhsIto(phi);
        phi(:,1) = y0;
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
        
        if(m>1)
            dbBridge = dbrownianBridge(deltaT, 0, xi);
        else
            dbBridge = zeros(order,1);
        end
        %Brownian bridge
        countRhsEvaluations = countRhsEvaluations + 1;
        
    end %choice of initialization
    
    
    %maximum number of sweeps
    for j=1:max_iter
        
        % SDC scheme for integration over all collocation points
        for p=2:points
            niter = 0;
            d_j = d(:,p-1);
            % deterministic part - 3fct evals
            rhs = rhs_eval(phi);
            rhs_integrate = (S(p-1,:)*rhs')';
            rhs_err = rhs_eval(phi(:,p-1)+d_j);
            rhs_phi_prev = rhs(:,p-1);
            countRhsEvaluations = countRhsEvaluations + 2;
            
            %stoch part
            stoch_rhs_integrate_b0 = zeros(order, 1);
            stoch_rhs_integrate_bm = zeros(order, 1);
            stoch_rhs_diff = zeros(order, 1);
            Jac_x = zeros(order, 1);
            for n=1:size(beta,2)
                %terms for stoch spectral integration
                stoch_rhs = stochRhs{n}(  phi());
                stoch_rhs_err = stochRhs{n}(  phi(:,p-1)+d_j );
                stoch_rhs_prev = stoch_rhs(:,p-1);
                countRhsEvaluations = countRhsEvaluations + 2;
                
                stoch_rhs_integrate_b0 = stoch_rhs_integrate_b0 + db0(n)*(S(p-1,:)*stoch_rhs')';
                if (m>1)
                    % loop over expansion terms: sum_1^m {xi_k * Q_k}
                    % include random variables xi as in the Karhunen-Loeve expansion
                    Q = zeros(1,size(quadMatK_c, 2));
                    for l=1:m-1
                        %Q = Q + cos(l*pi*t_currInt(1)/deltaT)*quadMatK_c(p-1,:, l)*xi(:,l) - sin(l*pi*t_currInt(1)/deltaT)*quadMatK_s(p-1,:, l)*xi(:,l);
                        Q = Q + cos(l*pi*t_currInt(1)/deltaT)*quadMatK_c(p-1,:, l)*xi(:,l);
                    end
                    Q= sqrt(2/deltaT)*Q;
                    stoch_rhs_integrate_bm = (Q*stoch_rhs')';
                    
                    dbBridge = dbrownianBridge(deltaT, t_currInt(p-1), xi);
                else
                    stoch_rhs_integrate_bm = zeros(order,1);
                end
                countRhsEvaluations = countRhsEvaluations + 1;
                
                stoch_rhs_diff = stoch_rhs_diff + (db0(n) + dbBridge(n))*(stoch_rhs_err - stoch_rhs_prev);
                Jac_x(n) = J{n}(phi(:,p-1));
                countRhsEvaluations = countRhsEvaluations + 1;
                
            end %end loop system's order
            
            %correction scheme
            F = - 1*((ones(length(p),1) - (t(:,p) - t(:,p-1)) * (db0(:) + dbBridge(:)) *Jac_x) *d(:,p)) + d(:,p-1)+(t(:,p) - t(:,p-1)) *(db0(:) + dbBridge(:))* ( rhs_err - rhs_phi_prev)...
                + rhs_integrate - phi(:,p) + phi(:, p-1) ...
                + stoch_rhs_integrate_b0 + stoch_rhs_integrate_bm;
            
            dF = ones(length(p),1) - (t(:,p) - t(:,p-1)) * (db0(:) + dbBridge(:)) *Jac_x;
            
            d_delta = linsolve(dF,F);
            d_j = d_j + d_delta;
            
            while abs(d_delta)>1e-15 && niter < 10
                
                %stoch part
                stoch_rhs_integrate_b0 = zeros(order, 1);
                stoch_rhs_integrate_bm = zeros(order, 1);
                stoch_rhs_diff = zeros(order, 1);
                Jac_x = zeros(order, 1);
                for n=1:size(beta,2)
                    %terms for stoch spectral integration
                    stoch_rhs = stochRhs{n}(  phi());
                    stoch_rhs_err = stochRhs{n}(  phi(:,p-1)+d_j );
                    stoch_rhs_prev = stoch_rhs(:,p-1);
                    countRhsEvaluations = countRhsEvaluations + 2;
                    
                    stoch_rhs_integrate_b0 = stoch_rhs_integrate_b0 + db0(n)*(S(p-1,:)*stoch_rhs')';
                    if (m>1)
                        stoch_rhs_integrate_bm = (Q*stoch_rhs')';
                    else
                        stoch_rhs_integrate_bm = zeros(order,1);
                    end
                    stoch_rhs_diff = stoch_rhs_diff + (db0(n) + dbBridge(n))*(stoch_rhs_err - stoch_rhs_prev);
                    countRhsEvaluations = countRhsEvaluations + 4;
                end %end loop system's order
                
                %correction scheme
                F = - 1*((ones(length(p),1) - (t(:,p) - t(:,p-1)) * (db0(:) + dbBridge(:)) *Jac_x) *d_j) + d(:,p-1)+(t(:,p) - t(:,p-1)) *(db0(:) + dbBridge(:))* ( rhs_err - rhs_phi_prev)...
                    + rhs_integrate - phi(:,p) + phi(:, p-1) ...
                    + stoch_rhs_integrate_b0 + stoch_rhs_integrate_bm;
                
                dF = ones(length(p),1) - (t(:,p) - t(:,p-1)) * (db0(:) + dbBridge(:)) *Jac_x;
                
                d_delta = linsolve(dF,F);
                d_j = d_j + d_delta;
                niter = niter+1;
            end %end while Newton scheme
            d(:,p)=d_j;
        end
        
        %solution update
        phi = phi + d;
        
        correctionNorm = norm(d);
        
        if (correctionNorm < tol)
            break;
        end
        
        %         if j==max_iter
        %             %correctionNorm
        %             fprintf('\n %d interval of %d: %g \n', i, intervals, correctionNorm)
        %         end
        
    end %iterations
    y0 = phi(:,end);
    sol(:,int_begin:int_end) = phi;
    d = zeros(order, points);
    t_begin = t_currInt(end);
    
end %intervals

end %function
