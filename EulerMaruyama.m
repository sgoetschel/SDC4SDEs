%% Euler Maruyama
% using the explicit Euler scheme
%by Lisa Fischer

function [ sol, countRhsEvaluations, time ] = EulerMaruyama( initial, deltaT, t_begin, t_end, rhs_eval, stoch_rhs, deltaW, order)

t = t_begin;
intervals = round((t_end-t_begin)/deltaT);
time = zeros(1,intervals);
time(1) = t_begin;

% initial value of ODE
phi= initial;
sol = zeros(order,intervals);
sol(:,1) = initial;

i=1;
countRhsEvaluations = 0;

while  (t+1e-9 <= t_end) 
    
    rhs_phi = rhs_eval(phi);
    countRhsEvaluations = countRhsEvaluations +1;
    
    stoch_part = zeros(order, 1);
    for n=1:order
        stochRhs = stoch_rhs{n}(phi);
        stoch_part = stoch_part + deltaW(n,i).*stochRhs;
        countRhsEvaluations = countRhsEvaluations +1;
    end
    
    %Euler-Maruyama scheme
    phi = phi + deltaT*rhs_phi   + stoch_part;
    
    t = t+deltaT;
    
    %save solutions and time points
    sol(:,i+1) = phi;
    time(i+1) = t;
    i= i+1;
    
end %intervals
end %function

