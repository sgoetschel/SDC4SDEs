function[rhs, stochRhs, J, RhsIto, exact] = problem1(lambda, beta)

rhs =@(x)(lambda-0.5*beta^2)*x;
stochRhs ={@(x) beta*x};
J = {@(x) beta};
RhsIto =@(x) lambda*x;

%exact sol of the Stratonovich SDE
exact =@exactExp;

end
