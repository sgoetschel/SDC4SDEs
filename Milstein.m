%% Milstein Method
%by Lisa Fischer
function [ sol, countRhsEvaluations, time ] = Milstein( initial, deltaT, t_begin, t_end, beta, rhs_eval, stoch_rhs, J, deltaW, order)

t = t_begin;
intervals = round((t_end-t_begin)/deltaT);
time = zeros(1,intervals+1);
time(1) = t_begin;

% initial value of ODE
phi= initial;
sol = zeros(order,intervals+1);
sol(:,1) = initial;

i=1;
countRhsEvaluations = 0;
dimBM = size(beta,2);


while  (t+1e-9 <= t_end)
    
    rhs_phi = rhs_eval(phi);
    countRhsEvaluations = countRhsEvaluations +1;
    
    stoch_part = [];
    %2-dim Milstein for commutative noise
    if (dimBM > 1)
        stochRhs = cell(1, dimBM);
        dstochRhs = cell(1, dimBM);
        for n=1:dimBM
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
    
    %Milstein scheme
    phi = phi + deltaT*rhs_phi...
        + stoch_part;
    t = t+deltaT;
    
    %save solutions and time points
    sol(:,i+1) = phi;
    time(i+1) = t;
    i= i+1;
    
end %intervals
end %function

