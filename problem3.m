function[rhs, stochRhs, J, RhsIto, exact] = problem3(lambda, beta)
%non-commutative noise
rhs =@(x) lambda*x.*(1-x)-0.5*beta*x;
stochRhs ={@(x) beta*x};
J = {@(x) beta};
RhsIto =@(x) lambda*x.*(1-x);

exact =@exactTP3;
end