%% Implicit Milstein scheme using damped Newton's method 
%by Lisa Fischer
function [ sol, countRhsEvaluations, time ] = ImplMilsteinDampNM( initial, deltaT, t_begin, t_end, beta, rhs_eval, stoch_rhs, J, deltaW, order)

t = t_begin;
intervals = round((t_end-t_begin)/deltaT);
time = zeros(1,intervals);
time(1) = t_begin;

% sol(:,i) value of ODE
phi= initial;
sol = zeros(order,intervals);
sol(:,1) = initial;

i=1;
countRhsEvaluations = 0;


while  (t+1e-9 <= t_end) %(t +deltaT <= t_end)
    
    rhs_phi = rhs_eval(phi);
    countRhsEvaluations = countRhsEvaluations +1;
    
    stoch_part = zeros(order, 1);
    %2-dim Milstein for commutative noise
    if (size(beta,2) > 1)
                for n=1:size(beta,2)
                    stochRhs{n} = stoch_rhs{n}(phi);
                    dstochRhs{n} = J{n}(phi);
                    countRhsEvaluations = countRhsEvaluations +2;
                end
        
        stoch_part = 0.5* ( (stochRhs{1}(1,:) * deltaW(1,i)^2 + stochRhs{2}(1,:) *deltaW(2,i)*deltaW(1,i))*dstochRhs{1}(:,1) + ...
                            (stochRhs{1}(2,:) * deltaW(1,i)^2 + stochRhs{2}(2,:)*deltaW(2,i)*deltaW(1,i))*dstochRhs{1}(:,2) + ...
                            (stochRhs{1}(1,:) * deltaW(1,i)*deltaW(2,i) + stochRhs{2}(1,:) *deltaW(2,i)^2)*dstochRhs{2}(:,1) +...
                            (stochRhs{1}(2,:) * deltaW(1,i)*deltaW(2,i) + stochRhs{2}(2,:) *deltaW(2,i)^2)*dstochRhs{2}(:,2) );
        
        stoch_part = stoch_part + [ stochRhs{1}, stochRhs{2} ]*deltaW(:,i);
    else %1-dim problem
        stochRhs = stoch_rhs{1}(phi);
        dstochRhs = J{1}(phi);
        stoch_part = deltaW(1,i).*stochRhs + 0.5* (stochRhs'*dstochRhs)'*(deltaW(1,i)^2-deltaT);
        countRhsEvaluations = countRhsEvaluations +2;
    end
      
    %implicit Milstein scheme using damped Newton's method    
        damping = 1;
        delta_phi = -(phi - rhs_phi*deltaT-sol(:,i)-stoch_part)/(1-(3*phi*phi-1)*deltaT);
        niter = 0;
        while abs(delta_phi)>1e-5 && niter < 10 
            phi_new = phi + damping*delta_phi;
            rhs_phi = rhs_eval(phi_new);
            countRhsEvaluations = countRhsEvaluations + 1;
            delta_phi_s = -(phi_new - rhs_phi*deltaT-sol(:,i)-stoch_part)/(1-(3*phi*phi-1)*deltaT);
            while abs(delta_phi_s)^2 > (1-damping/2)*abs(delta_phi)^2 && damping > 1e-6
                damping = 0.5*damping;
                phi_new = phi + damping*delta_phi;
                rhs_phi = rhs_eval(phi_new);
                countRhsEvaluations = countRhsEvaluations + 1;
                delta_phi_s = -(phi_new - rhs_phi*deltaT-sol(:,i)-stoch_part)/(1-(3*phi*phi-1)*deltaT);
            end
            damping = min(1,2*damping);
            phi = phi_new;
            delta_phi = -(phi - rhs_phi*deltaT-sol(:,i)-stoch_part)/(1-(3*phi*phi-1)*deltaT);
            niter = niter+1;
        end
    
    t = t+deltaT;
    
    %save solutions and time points
    sol(:,i+1) = phi;
    time(i+1) = t;
    i= i+1;
    
end %intervals
end %function

