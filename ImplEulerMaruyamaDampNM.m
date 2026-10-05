%% Implicit Euler Maruyama using damped Newton's method
%by Lisa Fischer

function [ sol, countRhsEvaluations, time ] = ImplEulerMaruyamaDampNM( initial, deltaT, t_begin, t_end, rhs_eval, stoch_rhs, deltaW, order)

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

while  (t+1e-9 <= t_end)
    
    rhs_phi = rhs_eval(phi);
    countRhsEvaluations = countRhsEvaluations +1;
    
    stoch_part = zeros(order, 1);
    for n=1:order
        stochRhs = stoch_rhs{n}(phi);
        stoch_part = stoch_part + deltaW(n,i).*stochRhs;
        countRhsEvaluations = countRhsEvaluations +1;
    end
    
    %implicit EM scheme using damped Newton's method
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

