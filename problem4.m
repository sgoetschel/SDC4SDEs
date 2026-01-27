function[rhs, stochRhs, J, RhsIto, exact] = problem4(A, B1, B2)

rhs =@(x) (A - 0.5*(B1*B1+B2*B2) )*x;
stochRhs = {@(x) B1*x, @(x) B2*x};
J = {@(x) B1 , @(x)B2};
RhsIto =@(x) A*x;

exact =@exactTP4;
end