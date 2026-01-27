function[rhs, stochRhs, J, RhsIto, exact] = problemOU(lambda, beta)

rhs =@(x)lambda*x;
stochRhs ={@(x) beta*ones(1,length(x))};
J = {@(x) 0};
RhsIto =@(x) lambda;

%exact sol of the Stratonovich SDE
exact = @exactOU;
end
